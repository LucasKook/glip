"""CLI: python -m glipcut <in.json> <out.json>. Options are read from the
"options" field of the input; errors are written to the output as "error".
With "task": "connected", the input holds a graph ("graph", in the start
format of the class) and queries; the output is the connection indicator
of every query (no solver needed)."""
import json
import sys
import traceback

from . import classes


def connected(case):
    G = classes.get(case.get('class', 'admg'))
    d = int(case['d'])
    val = G.start_values(d, case['graph'])
    return {'valid': not G.violated_cycle(d, val),
            'connected': [int(G.connected_vals(d, val, q['i'], q['j'], q['C']))
                          for q in case['queries']]}


def main(argv):
    if len(argv) != 3:
        sys.stderr.write('usage: python -m glipcut <in.json> <out.json>\n')
        return 2
    case = json.load(open(argv[1]))
    try:
        if case.get('task') == 'connected':
            res = connected(case)
        else:
            from .solver import solve
            res = solve(case)
        code = 0
    except Exception as ex:
        res = {'error': f'{type(ex).__name__}: {ex}', 'traceback': traceback.format_exc()}
        code = 1
    with open(argv[2], 'w') as f:
        json.dump(res, f)
    return code


if __name__ == '__main__':
    sys.exit(main(sys.argv))
