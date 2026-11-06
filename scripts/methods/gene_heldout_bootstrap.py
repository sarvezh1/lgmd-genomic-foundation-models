"""Reference gene-stratified paired locus bootstrap for pooled metrics."""
from __future__ import annotations
import argparse,csv,json,sys
from collections import defaultdict
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.dont_write_bytecode=True;sys.path.insert(0,str(ROOT/"scripts/validation"))
from metrics import auroc,average_precision,classification_metrics
def read_tsv(path):
 with Path(path).open(encoding="utf-8-sig",newline="") as h:return list(csv.DictReader(h,delimiter="\t"))
def write_tsv(path,rows):
 path.parent.mkdir(parents=True,exist_ok=True)
 with path.open("w",encoding="utf-8",newline="") as h:w=csv.DictWriter(h,fieldnames=list(rows[0]),delimiter="\t",lineterminator="\n");w.writeheader();w.writerows(rows)
def write_json(path,value):
 path.parent.mkdir(parents=True,exist_ok=True)
 with path.open("w",encoding="utf-8",newline="\n") as h:h.write(json.dumps(value,indent=2,sort_keys=True)+"\n")
def sample_indices(base,rng):
 import numpy as np
 sampled=[]
 for gene in sorted({r["GeneSymbol"] for r in base}):
  gene_rows=[(i,r) for i,r in enumerate(base) if r["GeneSymbol"]==gene];loci=sorted({r["locus_id"] for _,r in gene_rows});lookup=defaultdict(list)
  for i,r in gene_rows:lookup[r["locus_id"]].append(i)
  for locus in rng.choice(loci,size=len(loci),replace=True):sampled.extend(lookup[locus])
 return sampled
def calculate(records,indices):
 y=[int(records[i]["binary_label"]) for i in indices];p=[int(records[i]["prediction"]) for i in indices];s=[float(records[i]["score"]) for i in indices];c=classification_metrics(y,p)
 return c["mcc"],auroc(y,s),average_precision(y,s)
def run_bootstrap(rows,reps=1000,seed=1605):
 import numpy as np
 by=defaultdict(list)
 for r in rows:by[r["configuration_id"]].append(r)
 configs=sorted(by);base=by[configs[0]];rng=np.random.default_rng(seed);values={c:{"MCC":[],"AUROC":[],"AUPRC":[]} for c in configs}
 for _ in range(reps):
  indices=sample_indices(base,rng)
  for config in configs:
   observed=calculate(by[config],indices)
   for metric,value in zip(("MCC","AUROC","AUPRC"),observed):values[config][metric].append(value)
 intervals=[]
 model={r["configuration_id"]:r["model"] for r in rows}
 for config in configs:
  point=calculate(by[config],range(len(base)))
  for metric,estimate in zip(("MCC","AUROC","AUPRC"),point):
   intervals.append({"configuration_id":config,"model":model[config],"metric":metric,"point_estimate":estimate,"ci_lower":float(np.percentile(values[config][metric],2.5)),"ci_upper":float(np.percentile(values[config][metric],97.5)),"bootstrap_reps":reps,"seed":seed})
 return intervals
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument("--predictions",type=Path,default=ROOT/"data/gene_heldout_predictions.tsv");ap.add_argument("--archived-intervals",type=Path,default=ROOT/"data/gene_heldout_intervals.tsv");ap.add_argument("--output",type=Path,default=ROOT/"verification_outputs/gene_heldout_bootstrap/gene_heldout_bootstrap.tsv");ap.add_argument("--summary",type=Path,default=ROOT/"verification_outputs/gene_heldout_bootstrap/gene_heldout_bootstrap_summary.json");args=ap.parse_args()
 rows=read_tsv(args.predictions);observed=run_bootstrap(rows);arch={(r["model"],r["metric"]):r for r in read_tsv(args.archived_intervals)};maximum=0.0
 for r in observed:
  a=arch[(r["model"],r["metric"])];r["archived_ci_lower"]=a["ci_lower"];r["archived_ci_upper"]=a["ci_upper"];r["lower_abs_difference"]=abs(r["ci_lower"]-float(a["ci_lower"]));r["upper_abs_difference"]=abs(r["ci_upper"]-float(a["ci_upper"]));maximum=max(maximum,r["lower_abs_difference"],r["upper_abs_difference"]);r["status"]="PASS" if max(r["lower_abs_difference"],r["upper_abs_difference"])<=1e-12 else "FAIL"
 write_tsv(args.output,observed);summary={"implementation_status":"independently reconstructed reference method","design":"paired resampling of sorted genomic loci within sorted held-out genes","replicates":1000,"seed":1605,"interval_percentiles":[2.5,97.5],"configurations":8,"intervals":24,"maximum_absolute_difference":maximum,"all_archived_intervals_reproduced":all(r["status"]=="PASS" for r in observed)}
 write_json(args.summary,summary);print(json.dumps(summary,indent=2,sort_keys=True));return 0 if summary["all_archived_intervals_reproduced"] else 1
if __name__=="__main__":raise SystemExit(main())
