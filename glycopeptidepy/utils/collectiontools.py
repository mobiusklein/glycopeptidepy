from collections import defaultdict

try:
    from collections.abc import Mapping
except ImportError:
    from collections import Mapping

import itertools

try:
    range = xrange
except Exception:
    pass


def _identity(i):
    return i


def groupby(ungrouped_list, key_fn=_identity, transform_fn=_identity):
    groups = defaultdict(list)
    for item in ungrouped_list:
        key_value = key_fn(item)
        groups[key_value].append(transform_fn(item))
    return groups


def descending_combination_counter(counter):
    keys = counter.keys()
    count_ranges = map(lambda x: range(x + 1), counter.values())
    for combination in itertools.product(*count_ranges):
        yield dict(zip(keys, combination))


class decoratordict(dict):
    def __call__(self, key):
        def wrapper(f):
            self[key] = f
            return f
        return wrapper


class _AccumulatorBag(Mapping):
    def __init__(self, source=None):
        self.store = defaultdict(int)
        if source is not None:
            self.store.update(source)

    def __getitem__(self, i):
        return self.store[i]

    def __setitem__(self, i, v):
        self.store[i] = v

    def __add__(self, other):
        new = _AccumulatorBag(self)
        for key, value in other.items():
            new[key] += value
        return new

    def __iadd__(self, other):
        for key, value in other.items():
            self[key] += value
        return self

    def __iter__(self):
        return iter(self.store)

    def __len__(self):
        return len(self.store)

    def __contains__(self, key):
        return key in self.store

    def keys(self):
        return self.store.keys()

    def values(self):
        return self.store.values()

    def items(self):
        return self.store.items()

    def __repr__(self):
        return "_AccumulatorBag(%r)" % dict(self.store)


try:
    _descending_combination_counter = descending_combination_counter
    from glycopeptidepy._c.collectiontools import descending_combination_counter
except ImportError:
    has_c = False

try:
    from glycopeptidepy._c.count_table import CountTable

    del _AccumulatorBag

    class _AccumulatorBag(CountTable):
        pass
except ImportError:
    has_c = False


__all__ = [
    "descending_combination_counter", "_AccumulatorBag"
]
