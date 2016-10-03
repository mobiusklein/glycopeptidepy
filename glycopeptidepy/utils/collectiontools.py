from collections import defaultdict
import itertools

try:
    range = xrange
except Exception, e:
    pass

try:
    import cPickle as pickle
except:
    import pickle as pickle


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


def flatten(iterable):
    return tuple(itertools.chain.from_iterable(iterable))


class decoratordict(dict):
    def __call__(self, key):
        def wrapper(f):
            self[key] = f
            return f
        return wrapper
