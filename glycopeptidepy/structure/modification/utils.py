class ModificationIndex(dict):
    def __getitem__(self, key):
        try:
            return dict.__getitem__(self, key)
        except KeyError:
            return 0


class ModificationStringParseError(ValueError):
    pass


class ModificationNameResolutionError(KeyError):
    pass
