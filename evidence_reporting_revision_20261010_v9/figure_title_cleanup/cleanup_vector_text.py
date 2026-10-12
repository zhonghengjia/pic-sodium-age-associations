"""Approved, bounded removal of PDF text glyphs. No plot, estimate or model is recomputed."""
from pathlib import Path
import hashlib,json,shutil,math
import pymupdf as fitz
import numpy as np
from PIL import Image
R=Path(__file__).resolve().parent;P=R.parent.parent
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
rows=json.loads((R/'inventory.json').read_text(encoding='utf8'))
source_R={
'Figure1':'figures_tables_mainline_20261009_v1/code/render_display.R',
'Figure2':'analysis_daily_mean_20261009_v1/figure_association_roundticks_20261009_v10/code/render_round_ticks.R',
'Figure3':'side_figure3_dynamic_composite_20261010_v1/code/Figure3_dynamic_composite.R',
'Figure4':'side_day1_day3_ROC_age_20261010_v2/code/render_two_models.R',
'Figure5':'KM_main_annotation_removed_20261009_v7/code/draw_main.R',
'Supplementary_Figure1':'figure2_language_revision_20261010_v4/code/render_supplement_figures.R',
'Supplementary_Figure2':'supplementary_figure2_caption_20261010_v6/code/render_supplementary_figure2.R',
'Supplementary_Figure3':'figure2_language_revision_20261010_v4/code/render_supplement_figures.R',
'Supplementary_Figure4':'sodium_day3_CIF_color_revision_20261009_v2/code/render.R'}
rules={
'Figure1':['Cohort selection'],
'Figure2':['Age interaction Holm P='],
'Figure3':['A  Observed three-day sodium changes','B  Adjusted associations of sodium changes','Unadjusted group means','OR per 5 mmol/L greater change','Survived n='],
'Figure4':['Day 1 versus day 3 sodium: ROC by age','Clinical models evaluated','n = '],
'Figure5':['Day-3 sodium and subsequent in-hospital survival across ages','Unadjusted Kaplan-Meier estimates','N = ','Three-group log-rank'],
'Supplementary_Figure1':[], 'Supplementary_Figure2':[],
'Supplementary_Figure3':['Intercept '],
'Supplementary_Figure4':['Day 3 sodium and subsequent in-hospital mortality','n = ']}
for sub in ['original','clean','comparison','source_code']:(R/sub).mkdir(exist_ok=True)
receipts=[];moves=[];source_manifest=[]
for row in rows:
    key=row['id'];src=P/row['pdf_matches'][0];orig=R/'original'/f'{key}.pdf';shutil.copyfile(src,orig)
    code=P/source_R[key];shutil.copyfile(code,R/'source_code'/f'{key}_original.R')
    source_manifest.append(dict(figure=key,original_pdf=str(src),pdf_sha256=sha(src),original_embedded_png_sha256=row['sha256'],R_entry=str(code),R_entry_sha256=sha(code)))
    doc=fitz.open(orig);page=doc[0];rects=[];deleted=[];retained=[]
    original_rect=list(page.rect);original_drawings=page.get_drawings();original_images=page.get_images(full=True)
    raw=page.get_text('rawdict')
    for block in raw['blocks']:
        for line in block.get('lines',[]):
            for span in line['spans']:
                chars=span['chars'];text=''.join(c['c'] for c in chars)
                if not any(text.startswith(pat) for pat in rules[key]):continue
                start=0
                if key=='Figure3' and text.startswith(('A  Observed','B  Adjusted')):
                    # Keep the existing panel letter and its location. Only following descriptive words leave.
                    start=3
                targets=chars[start:];r=fitz.Rect(targets[0]['bbox'])
                for c in targets[1:]:r|=fitz.Rect(c['bbox'])
                # Text-only deletion: no white fill, no replacement text, no line/image redaction.
                rects.append(list(r));deleted.append(text[start:])
                moves.append(dict(figure=key,exact_original=text,removed_text=text[start:],retained_panel_prefix=text[:start],bbox_pt=list(r),destination='existing figure legend/caption or existing population table by reference'))
    expected={'Figure1':1,'Figure2':3,'Figure3':10,'Figure4':8,'Figure5':14,'Supplementary_Figure1':0,'Supplementary_Figure2':0,'Supplementary_Figure3':2,'Supplementary_Figure4':7}[key]
    assert len(deleted)==expected,(key,len(deleted),expected)
    clean=R/'clean'/f'{key}.pdf'
    if deleted:
        # Edit original literal text string operands only. This preserves every
        # original graphics/text coordinate byte; generic PDF redaction rewrites
        # retained text positioning and is deliberately not used.
        original_stream=page.read_contents();stream=original_stream
        def escape_literal(s):return s.encode('ascii').replace(b'\\',b'\\\\').replace(b'(',b'\\(').replace(b')',b'\\)')
        for s in sorted(set(deleted),key=len,reverse=True):
            old=escape_literal(s);count=stream.count(old)
            assert count==deleted.count(s),(key,'literal object not uniquely bound',s,count,deleted.count(s))
            stream=stream.replace(old,b'')
        streams=page.get_contents();assert len(streams)==1,(key,'multiple content streams')
        doc.update_stream(streams[0],stream)
        doc.save(clean,garbage=4,deflate=True)
    else:shutil.copyfile(orig,clean)
    doc.close()
    with fitz.open(orig) as old,fitz.open(clean) as new:
        before=old[0];after=new[0]
        assert list(after.rect)==original_rect,key+' canvas changed'
        # Display-list sequence indices change when a text draw becomes empty;
        # every scientific drawing property and its order must remain exact.
        drawing_geometry=lambda pg:[{k:v for k,v in a.items() if k!='seqno'} for a in pg.get_drawings()]
        assert drawing_geometry(before)==drawing_geometry(after),key+' vector drawing changed'
        assert before.get_images(full=True)==after.get_images(full=True),key+' embedded image changed'
        # Every retained glyph must retain character, exact bbox, origin, font, size and color.
        def glyphs(pg):
            out=[]
            for block in pg.get_text('rawdict')['blocks']:
                for line in block.get('lines',[]):
                    for sp in line['spans']:
                        for ch in sp['chars']:out.append((ch['c'],tuple(ch['bbox']),tuple(ch['origin']),sp['font'],sp['size'],sp['color']))
            return out
        oldglyphs=glyphs(before);newglyphs=glyphs(after)
        # Redaction may omit whitespace at removed regions; multiset retained nonspace glyphs must be exact.
        from collections import Counter
        nc=Counter(g for g in newglyphs if not g[0].isspace());oc=Counter(g for g in oldglyphs if not g[0].isspace())
        assert not (nc-oc),key+' new or shifted glyph'
        removed=oc-nc
        for g in removed:
            assert any(fitz.Rect(g[1]).intersects(fitz.Rect(r)) for r in rects),key+' glyph removed outside approved rect'
        # Compare equal-renderer PNGs at exactly the embedded original PNG dimensions.
        with Image.open(R/'original'/f'{key}.png') as im:w,h=im.size
        matrix=fitz.Matrix(w/before.rect.width,h/before.rect.height)
        p0=before.get_pixmap(matrix=matrix,alpha=False);p1=after.get_pixmap(matrix=matrix,alpha=False)
        assert (p0.width,p0.height)==(p1.width,p1.height)
        a=np.frombuffer(p0.samples,dtype=np.uint8).reshape(p0.height,p0.width,p0.n)
        b=np.frombuffer(p1.samples,dtype=np.uint8).reshape(p1.height,p1.width,p1.n)
        diff=np.any(a!=b,axis=2);mask=np.zeros(diff.shape,dtype=bool)
        for r in rects:
            rr=fitz.Rect(r)*matrix
            x0=max(0,math.floor(rr.x0)-2);y0=max(0,math.floor(rr.y0)-2)
            x1=min(p0.width,math.ceil(rr.x1)+2);y1=min(p0.height,math.ceil(rr.y1)+2)
            mask[y0:y1,x0:x1]=True
        outside=int(np.count_nonzero(diff&~mask));assert outside==0,(key,'pixel change outside approved text',outside)
        p0.save(R/'comparison'/f'{key}_before_pdf_render.png')
        p1.save(R/'comparison'/f'{key}_after_pdf_render.png')
        if deleted:p1.save(R/'clean'/f'{key}.png')
        else:shutil.copyfile(R/'original'/f'{key}.png',R/'clean'/f'{key}.png')
        text=after.get_text()
        for s in deleted:assert s not in text,(key,'removed text still present',s)
        if key=='Figure3':assert 'A' in text and 'B' in text
        receipts.append(dict(figure=key,removed_lines=len(deleted),unchanged_original_bytes=not deleted,canvas_pt=original_rect,
            original_docx=row['docx']+'.docx',original_docx_media=row['media'],original_pdf_source=str(src),new_png=str(R/'clean'/f'{key}.png'),new_pdf=str(clean),
            canvas_unchanged=True,vector_drawings_exact=True,embedded_images_exact=True,retained_nonspace_glyphs_exact=True,
            removed_glyph_count=sum(removed.values()),changed_pixels=int(np.count_nonzero(diff)),pixels_changed_outside_approved_text_boxes=outside,
            comparison_renderer='PyMuPDF identical before/after Matrix at original embedded PNG dimensions',
            original_png_render_backend='original R ragg; archived byte-for-byte',final_png_backend='clean PDF rasterized with PyMuPDF (approved vector deletion route)' if deleted else 'unchanged original R ragg PNG',
            warning='Comparison before PNG is the original PDF rendered by the same engine, not the original ragg PNG. Cairo/ragg antialiasing differs.',
            approved_text_rectangles_pt=rects,original_pdf_sha256=sha(orig),clean_pdf_sha256=sha(clean),clean_png_sha256=sha(R/'clean'/f'{key}.png'),
            models_refitted=False,data_read_or_recomputed=False,no_white_fill=True,no_raster_cover=True))
        print(key,len(deleted),'removed lines; vector drawings exact; outside-text pixel diff',outside)
(R/'source_manifest.json').write_text(json.dumps(source_manifest,ensure_ascii=False,indent=2),encoding='utf8')
(R/'text_relocation_manifest.json').write_text(json.dumps(moves,ensure_ascii=False,indent=2),encoding='utf8')
(R/'verification_receipt.json').write_text(json.dumps(receipts,ensure_ascii=False,indent=2),encoding='utf8')
md=['# Approved text relocation list','', 'Original data, estimates, curves, drawing objects, axes and retained glyph locations are unchanged. Only registered text was removed without white fill. Figures S1 and S2 remain byte-for-byte unchanged.','']
for row in rows:
    key=row['id'];md+=['## '+key,'']
    mm=[m for m in moves if m['figure']==key]
    if not mm:md+=['No change: no standalone figure title or statistical/method subtitle.','']
    for m in mm:md+=['- Exact original: `'+m['exact_original']+'`'+(' (panel letter retained)' if m['retained_panel_prefix'] else ''),'']
(R/'TEXT_RELOCATION.md').write_text('\n'.join(md),encoding='utf8')
