import re
from numpy import int32, int64

def parse_index(string, shift=0):
    if not string:
        return []
    assert re.match(r'(\d+(\-\d+)?)(,\d+(\-\d+)?)*', string)
    segs = string.strip().split(',')
    indices = []
    for seg in segs:
        m = re.match(r'(\d+)\-(\d+)', seg)
        if m:
            indices += list(range(int(m.group(1)), int(m.group(2)) + 1))
        else:
            indices.append(int(seg))
    return [i+shift for i in indices]

def rev_parse_index(l):
    templ = list(l)
    for ele in templ:
        assert type(ele) in (int, int32, int64)
    templ.sort()
    start = templ.pop(0)
    end = start
    s = ''
    while templ:
        ele = templ.pop(0)
        if ele == end:
            continue
        if ele == end + 1:
            end = ele
        else:
            if start == end:
                s = s + '%d,' % start
            else:
                s = s + '%d-%d,' % (start, end)
            start = ele
            end = start
    if start == end:
        s = s + '%d' % start
    else:
        s = s + '%d-%d' % (start, end)
    return s


