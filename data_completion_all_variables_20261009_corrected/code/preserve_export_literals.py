"""Preserve every nonmissing source literal without touching any imputed cell."""
from pathlib import Path
import csv,json,hashlib
O=Path(__file__).resolve().parents[1]
def sha(p):
    with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read(p):
    with p.open(encoding='utf-8-sig',newline='') as f:
        r=csv.DictReader(f);return r.fieldnames,list(r)
source=O/'data/restricted/clean_observed_patient.csv';fields,raw=read(source)
receipt=json.loads((O/'results/MI/execution_receipt.json').read_text(encoding='utf-8'))
assert receipt['source_hashes']['input']==sha(source) and receipt['source_hashes']['export_code']==sha(__file__)
records=[]
for j in range(1,51):
    p=O/f'results/restricted_MI/completed_{j:02}.csv';cols,rows=read(p)
    assert cols==fields and len(rows)==len(raw)==3459
    missing_before=[]
    for a,b in zip(raw,rows):
        assert a['id']==b['id']
        for c in cols:
            if a[c]!='':b[c]=a[c]
            else:missing_before.append((a['id'],c,b[c]))
    before_hash=hashlib.sha256(json.dumps(missing_before,ensure_ascii=False).encode()).hexdigest()
    with p.open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=cols);w.writeheader();w.writerows(rows)
    _,check=read(p);after=[]
    for a,b in zip(raw,check):
        for c in cols:
            if a[c]!='':assert a[c]==b[c]
            else:after.append((a['id'],c,b[c]))
    assert missing_before==after
    records.append({'imputation':j,'sha256':sha(p),'nonmissing_literals_exact':True,'original_missing_export_block_unchanged':True,'missing_block_sha256':before_hash})
(O/'results/MI/export_literal_validation.json').write_text(json.dumps({'source_sha256':sha(source),'script_sha256':sha(__file__),'M':50,'checks':records},ensure_ascii=False,indent=2),encoding='utf-8')
print('All 50 exported tables preserve original nonmissing literals; imputed blocks unchanged.')
