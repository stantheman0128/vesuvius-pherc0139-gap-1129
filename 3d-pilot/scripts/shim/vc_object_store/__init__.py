"""Stand-in for the compiled vc_object_store package. infer_rowcol_triplet_wraps.py imports it at module level
but uses it only in --store mode, which this smoke test does not use. Any use raises."""


class ObjectStore:
    @classmethod
    def open(cls, *a, **k):
        raise RuntimeError("shim vc_object_store: --store mode is not available here")


class _Surfaces:
    def __getattr__(self, name):
        raise RuntimeError(f"shim vc_object_store.surfaces.{name}: not available here")


surfaces = _Surfaces()
