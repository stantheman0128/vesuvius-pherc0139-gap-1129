"""Small metadata: zarr .zarray (levels 0-3), .zattrs, metadata.json, transform.json, umbilicus,
meta.json of the on-1.129 tifxyz copies of w034/w035/w036. Into downloads/meta/."""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from netlog import fetch, ROOT  # noqa: E402

BUCKET = "https://vesuvius-challenge-open-data.s3.amazonaws.com"
VOL = f"{BUCKET}/PHerc0139/volumes/20260413113053-1.129um-0.2m-59keV-masked.zarr"
SEG = {
    "w034": "20260303000000-w034_2026030317",
    "w035": "20260317000000-w035_2026031718",
    "w036": "20260316000000-w036_2026031607",
}
OUT = os.path.join(ROOT, "downloads", "meta")


def get(url, name):
    dst = os.path.join(OUT, name)
    if os.path.exists(dst):
        return
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    try:
        data = fetch(url, tag="meta")
    except Exception as exc:
        print("MISSING", url, exc)
        return
    open(dst, "wb").write(data)
    print(len(data), name)


for lv in range(4):
    get(f"{VOL}/{lv}/.zarray", f"zarr/{lv}.zarray")
get(f"{VOL}/.zattrs", "zarr/.zattrs")
get(f"{VOL}/metadata.json", "zarr/metadata.json")
get(f"{VOL}/meta.json", "zarr/meta.json")
get(f"{VOL}/transform.json", "zarr/transform.json")
get(f"{BUCKET}/PHerc0139/representations/umbilicus/20260102150214-umbilicus-20260826112529.json", "umbilicus.json")
for w, sid in SEG.items():
    ts = sid.split("-")[0]
    get(f"{BUCKET}/PHerc0139/segments/{sid}/mesh/{ts}-on-20260413113053-1.129um.tifxyz/meta.json", f"{w}-1.129-meta.json")
