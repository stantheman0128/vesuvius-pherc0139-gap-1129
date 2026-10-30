# Where the scripts read inputs and write outputs. Imported first by every script.
#   PHERC0139_DATA  download cache (catalog, umbilicus, meshes, CT chunks). Default: data/ at the repo root.
#   PHERC0139_OUT   generated outputs (.npz caches, PNGs, *_layers.json). Default: out/ at the repo root.
# The committed copies of the outputs are in results/ and png/.
import os, site, sys

# python -I leaves out the user site directory; packages installed with `pip install --user` live there.
sys.path.append(site.getusersitepackages())

REPO = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.abspath(os.environ.get("PHERC0139_DATA") or os.path.join(REPO, "data"))
OUT_ROOT = os.path.abspath(os.environ.get("PHERC0139_OUT") or os.path.join(REPO, "out"))
OUT_PILOT = os.path.join(OUT_ROOT, "pilot")
OUT_PLACES = os.path.join(OUT_ROOT, "places")
for _d in (DATA, OUT_PILOT, OUT_PLACES):
    os.makedirs(_d, exist_ok=True)

BUCKET = "https://vesuvius-challenge-open-data.s3.amazonaws.com"
# segment folder ids of the five windings (PHerc0139/segments/<id>/mesh/...)
SEGMENT_IDS = {
    "w033": "20260219000000-w033_2026021913",
    "w034": "20260303000000-w034_2026030317",
    "w035": "20260317000000-w035_2026031718",
    "w036": "20260316000000-w036_2026031607",
    "w037": "20260310000000-w037_2026031015",
}


def mesh_url(w, frame):
    """URL prefix of the tifxyz folder of mesh w; frame is '2.399' or '9.362'."""
    sid = SEGMENT_IDS[w]
    vol = {"2.399": "20260102150214-2.399um", "9.362": "20250728140407-9.362um"}[frame]
    return f"{BUCKET}/PHerc0139/segments/{sid}/mesh/{sid.split('-')[0]}-on-{vol}.tifxyz"
