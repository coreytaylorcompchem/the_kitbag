class BaseBackend:
    def __init__(self, **kwargs):
        self.pose = None
        self.cache = {}

    def load_pdb(self, pdb_path):
        raise NotImplementedError("Subclasses must implement `load_pdb` method.")

    def save_pose(self, path):
        if self.pose:
            self.pose.dump_pdb(path)
        else:
            raise RuntimeError("No pose is loaded in backend.")
