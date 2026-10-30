# Download the small inputs into $PHERC0139_DATA (default data/). Files already there are skipped.
#   python -I fetch_inputs.py              catalog, umbilicus, 9.362 um copies of w033-w037 (about 8 MB)
#   python -I fetch_inputs.py --coverage   also what coverage/ reads (volume metadata, 2.399 um mesh meta.json,
#                                          two level-5 chunks, the level-3 key listing, the #1381 transform)
# CT chunks and the 2.399 um mesh rows are fetched later by the pilot/ and places/ scripts, only where needed.
import os, sys, re, urllib.request, urllib.parse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from repo_paths import DATA, BUCKET, SEGMENT_IDS, mesh_url

VOL11 = f"{BUCKET}/PHerc0139/volumes/20260413113053-1.129um-0.2m-59keV-masked.zarr"
NB = [0]


def get(url, dest):
    dest = os.path.join(DATA, dest)
    if os.path.exists(dest):
        return
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    data = urllib.request.urlopen(url, timeout=120).read()   # the catalog is stored gzip-encoded; kept as is
    open(dest, "wb").write(data)
    NB[0] += len(data)
    print(f"  {len(data):>10,d}  {dest}")


def list_keys(prefix, dest):
    dest = os.path.join(DATA, dest)
    if os.path.exists(dest):
        return
    lines, tok = [], None
    while True:
        q = {"list-type": "2", "prefix": prefix}
        if tok:
            q["continuation-token"] = tok
        d = urllib.request.urlopen(f"{BUCKET}/?{urllib.parse.urlencode(q)}", timeout=120).read()
        NB[0] += len(d)
        d = d.decode("utf-8")
        lines += [f"{k} {s}" for k, s in re.findall(r"<Key>([^<]*)</Key>.*?<Size>([^<]*)</Size>", d)]
        m = re.search(r"<NextContinuationToken>([^<]*)</NextContinuationToken>", d)
        if not m:
            break
        tok = m.group(1)
    open(dest, "w").write("\n".join(lines) + "\n")
    print(f"  {len(lines):>10,d} keys  {dest}")


print(f"data dir: {DATA}")
get(f"{BUCKET}/metadata.json", "catalog-metadata.json")
get(f"{BUCKET}/PHerc0139/representations/umbilicus/20260102150214-umbilicus-20260826112529.json", "umbilicus.json")
for w in SEGMENT_IDS:
    for f in ("x.tif", "y.tif", "z.tif", "meta.json"):
        get(f"{mesh_url(w, '9.362')}/{f}", f"seg/{w}-9.362/{f}")

if "--coverage" in sys.argv:
    for w in SEGMENT_IDS:
        get(f"{mesh_url(w, '2.399')}/meta.json", f"seg/{w}-on-20260102150214-2.399um-meta.json")
    v = "20260413113053-1.129um-0.2m-59keV-masked"
    get(f"{VOL11}/0/.zarray", f"{v}/0_.zarray")
    get(f"{VOL11}/metadata.json", f"{v}/metadata.json")
    get(f"{VOL11}/transform.json", f"{v}/transform.json")
    for k in ("1/5/5", "3/4/6"):
        get(f"{VOL11}/5/{k}", f"L5/{k.replace('/', '_')}")
    get("https://raw.githubusercontent.com/7jycwjmbfn-eng/pherc0139-physical-audit/"
        "3f2bc75e63fbdabb3926e4d550b3a5851a231dfb/20260413113053-to-20250728140407.json",
        "pr1381-20260413113053-to-20250728140407.json")
    list_keys("PHerc0139/volumes/20260413113053-1.129um-0.2m-59keV-masked.zarr/3/", "keys-L3.txt")

print(f"downloaded {NB[0]:,d} bytes this run")
