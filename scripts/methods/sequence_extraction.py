"""Portable GRCh38 sequence extraction using a FASTA index.

This strict implementation follows the study accession, coordinate,
orientation, and alternate-substitution method. It fails on invalid
chromosomes, boundaries, or central-reference mismatches. Coordinates are
one-based; output windows are forward-strand genomic sequence.
"""
from __future__ import annotations
import argparse, csv, hashlib, json
from pathlib import Path

ACCESSIONS = {
 "1":"NC_000001.11", "2":"NC_000002.12", "3":"NC_000003.12",
 "4":"NC_000004.12", "5":"NC_000005.10", "6":"NC_000006.12",
 "7":"NC_000007.14", "8":"NC_000008.11", "9":"NC_000009.12",
 "10":"NC_000010.11", "11":"NC_000011.10", "12":"NC_000012.12",
 "13":"NC_000013.11", "14":"NC_000014.9", "15":"NC_000015.10",
 "16":"NC_000016.10", "17":"NC_000017.11", "18":"NC_000018.10",
 "19":"NC_000019.10", "20":"NC_000020.11", "21":"NC_000021.9",
 "22":"NC_000022.11", "X":"NC_000023.11", "Y":"NC_000024.10",
 "MT":"NC_012920.1", "M":"NC_012920.1",
}

def read_tsv(path):
 with Path(path).open(encoding="utf-8-sig",newline="") as h:return list(csv.DictReader(h,delimiter="\t"))
def sha256(path):
 d=hashlib.sha256()
 with Path(path).open("rb") as h:
  for block in iter(lambda:h.read(1048576),b""):d.update(block)
 return d.hexdigest()
def load_fai(path):
 index={}
 for line in Path(path).read_text(encoding="ascii").splitlines():
  name,length,offset,bases,width=line.split("\t")[:5];index[name]=(int(length),int(offset),int(bases),int(width))
 return index
def fetch(handle,record,start,end):
 length,offset,line_bases,line_width=record
 if start<0 or end>length or end<=start:raise ValueError(f"out-of-bounds interval {start}:{end} for contig length {length}")
 result=bytearray();pos=start
 while pos<end:
  within=pos%line_bases;take=min(end-pos,line_bases-within);byte=offset+(pos//line_bases)*line_width+within
  handle.seek(byte);result.extend(handle.read(take));pos+=take
 return result.decode("ascii").upper()
def construct(handle,index,chromosome,position,ref,alt,context):
 if context<=0 or context%2!=1:raise ValueError("context must be a positive odd integer")
 chrom=chromosome.removeprefix("chr");accession=ACCESSIONS.get(chrom)
 if not accession or accession not in index:raise ValueError(f"no primary-assembly accession for chromosome {chromosome}")
 half=context//2;center0=int(position)-1;sequence=fetch(handle,index[accession],center0-half,center0+half+1)
 if sequence[half]!=ref.upper():raise ValueError(f"reference mismatch at {chromosome}:{position}: FASTA={sequence[half]} table={ref}")
 alternate=sequence[:half]+alt.upper()+sequence[half+1:]
 return sequence,alternate,accession
def validate_context(fasta,fai,archived):
 rows=read_tsv(archived); context=len(rows[0]["reference_sequence"]); index=load_fai(fai); mismatches=[]; ambiguous=0
 with Path(fasta).open("rb") as handle:
  for row in rows:
   try:
    ref,alt,_=construct(handle,index,row["Chromosome"],row["Position"],row["REF"],row["ALT"],context)
   except Exception as error:
    raise RuntimeError(
     f'failed VariationID={row["VariationID"]} chromosome={row["Chromosome"]} '
     f'position={row["Position"]} context={context}'
    ) from error
   ambiguous+=sum(base not in "ACGT" for base in ref)
   if ref!=row["reference_sequence"] or alt!=row["alternate_sequence"]:mismatches.append(row["VariationID"])
 return {"context_bp":context,"records":len(rows),"exact_reference_matches":len(rows)-len(mismatches),"exact_alternate_matches":len(rows)-len(mismatches),"mismatch_count":len(mismatches),"first_mismatch_ids":mismatches[:10],"ambiguous_reference_bases":ambiguous,"archived_sha256":sha256(archived)}
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument("--fasta",type=Path,required=True);ap.add_argument("--fai",type=Path,required=True);ap.add_argument("--archived",type=Path,nargs="+",required=True);ap.add_argument("--report",type=Path,required=True);args=ap.parse_args()
 results=[validate_context(args.fasta,args.fai,p) for p in args.archived]
 report={"implementation_status":"strict portable study implementation","reference_accession":"GCF_000001405.40","assembly":"GRCh38.p14","coordinate_system":"one-based input converted to zero-based half-open FASTA intervals","orientation":"forward genomic strand","boundary_policy":"fail on out-of-bounds interval","contexts":results,"all_exact":all(r["mismatch_count"]==0 for r in results),"fasta_fai_sha256":sha256(args.fai)}
 args.report.parent.mkdir(parents=True,exist_ok=True);args.report.write_text(json.dumps(report,indent=2,sort_keys=True)+"\n",encoding="utf-8");print(json.dumps(report,indent=2,sort_keys=True));return 0 if report["all_exact"] else 1
if __name__=="__main__":raise SystemExit(main())
