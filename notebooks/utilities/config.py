class Config(dict):
    def __setattr__(self, name, value):
        self[name]  = value