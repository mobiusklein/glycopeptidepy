import unittest

from glycopeptidepy.utils import collectiontools


class TestCollectionTools(unittest.TestCase):
    def test_groupby(self):
        in_data = ["SPAM", "SPAM", "MAPS"]
        groups = collectiontools.groupby(in_data)
        assert len(groups) == 2
        assert sum({k: len(v) for k, v in groups.items()}.values()) == len(in_data)
        # test key_fn changing aggregation
        groups = collectiontools.groupby(in_data, key_fn=frozenset)
        assert len(groups) == 1
        assert sum({k: len(v) for k, v in groups.items()}.values()) == len(in_data)

    def test_accumulator_bag(self):
        bag = collectiontools._AccumulatorBag({"a": 4, "b": 3, (1, 2, 3): 2})
        abag = collectiontools._AccumulatorBag({"a": 8, "b": 6, (1, 2, 3): 4})
        bbag = bag + bag
        assert abag == bbag
        assert len(bag) == 3
        assert "a" in bag
        assert (1, 2, 3) in bag
        assert list(bag) == list(bag.keys())
        assert sum(bag.values()) == 4 + 3 + 2


if __name__ == '__main__':
    unittest.main()
