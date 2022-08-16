import sys
import json
from functools import partial

def add_names(mod_dict, alt_names=None):
    if alt_names is None:
        alt_names = []
    for name in alt_names:
        if name not in mod_dict['alt_names']:
            mod_dict['alt_names'].append(name)

def remove_names(mod_dict, alt_names=None):
    if alt_names is None:
        alt_names = []
    for name in alt_names:
        try:
            mod_dict['alt_names'].remove(name)
        except ValueError as err:
            print(err, name, mod_dict['alt_names'])


patches = [
    {
        "match_key": "full_name",
        "match_value": "Dehydration",
        "changes": [partial(add_names, alt_names=['Dehydrated']), partial(remove_names, alt_names=["Phospho+PL"])]
    }
]


def process(modifications):
    for patch in patches:
        for mod_dict in modifications:
            if mod_dict[patch['match_key']] == patch['match_value']:
                for change in patch['changes']:
                    change(mod_dict)

def main():
    infile = sys.argv[1]
    try:
        outfile = sys.argv[2]
    except IndexError:
        outfile = infile
    with open(infile, 'rt') as fh:
        modification_dicts = json.load(fh)
    process(modification_dicts)
    with open(outfile, 'wt') as fh:
        json.dump(modification_dicts, fh, sort_keys=True, indent=4)

if __name__ == "__main__":
    main()
