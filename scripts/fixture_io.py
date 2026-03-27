import gzip
import json
import os


def resolve_fixture_path(path: str) -> str:
    if os.path.exists(path):
        return path

    gz_path = f"{path}.gz"
    if os.path.exists(gz_path):
        return gz_path

    raise FileNotFoundError(f"fixture not found: {path} or {gz_path}")


def load_json_fixture(path: str):
    resolved = resolve_fixture_path(path)
    if resolved.endswith(".gz"):
        with gzip.open(resolved, "rt") as handle:
            return json.load(handle)

    with open(resolved) as handle:
        return json.load(handle)


def dump_json_gz(path: str, payload) -> None:
    with gzip.open(path, "wt") as handle:
        json.dump(payload, handle)