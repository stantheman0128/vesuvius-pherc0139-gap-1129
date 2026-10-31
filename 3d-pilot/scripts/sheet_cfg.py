"""Print the sheet tracer checkpoint's data/trace/model config (loaded with weights_only=True) and the files the
config refers to. Usage (WSL venv): python -I sheet_cfg.py <last.pt> <out_json>"""
import json
import os
import sys

os.environ["OMP_NUM_THREADS"] = "1"
import torch  # noqa: E402

torch.set_num_threads(1)
ck = torch.load(sys.argv[1], map_location="cpu", weights_only=True)
print("top keys:", list(ck.keys()))
cfg = ck["config"]
for sec in ("data", "trace", "model"):
    print(f"[{sec}]")
    for k, v in (cfg.get(sec) or {}).items():
        s = json.dumps(v, default=str)
        print(f"  {k} = {s[:160]}")
json.dump(cfg, open(sys.argv[2], "wt"), indent=1, default=str)
