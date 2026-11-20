#!/bin/sh
# igv-reports embeds igv.js into its pages and would download it in every job; fetch the
# version its template names once per environment instead, so the report jobs need no network
set -eu
python - << '_EOF'
import hashlib
import os
import urllib.request

url = "https://cdn.jsdelivr.net/npm/igv@3.5.2/dist/igv.min.js"
sha256 = "0efd638a0997aa90791ce6c83a8b33912d4bc06aed5e28741fc74f32a20998d6"
asset = os.path.join(os.environ["CONDA_PREFIX"], "share", "igv-reports", "igv.min.js")
os.makedirs(os.path.dirname(asset), exist_ok=True)
with urllib.request.urlopen(url) as response:
    data = response.read()
if hashlib.sha256(data).hexdigest() != sha256:
    raise SystemExit(f"{url}: checksum mismatch")
with open(asset, "wb") as handle:
    handle.write(data)
_EOF
