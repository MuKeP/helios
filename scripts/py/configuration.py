
import sys
from itertools import product, permutations

####################################################################################################

def sort_order(mult):
    label, occ, vac = mult

    # most significant amplitude kind (...,r,t)
    k = ord(label)-ord('a')

    # next - order of excitation (8,...,3,2,1), greater -> first
    f = -len(occ)

    # next if same amplitudes (e.g. t2 and t2) found, sort by occupied
    o1 = occ[0]

    # and second occupied, if exist
    o2 = 0
    if len(occ) > 1:
        o2 = occ[1]

    # and third occupied, if exist
    o3 = 0
    if len(occ) > 2:
        o3 = occ[2]

    return 10000*k+1000*f+100*o1+10*o2+1*o3

####################################################################################################

def configuration_split(occ, vac, separation, labels):
    component = []

    iterator = iter(zip(occ, vac))
    for label, fold in zip(labels, separation):

        # pop occupied and vacant simultaneously
        occ_, vac_ = [], []
        for _ in range(fold):
            o, v = next(iterator)
            occ_.append(o)
            vac_.append(v)

        # since we use spin-orbital basis, it is allowed to resort indices without constrains
        component.append((label, tuple(sorted(occ_)), tuple(sorted(vac_))))
    return frozenset(component)

####################################################################################################

def generate_components(kind, expected={'t', 'r'}, verbose=True):

    amplitudes = kind.split('*')

    labels = [a[0] for a in amplitudes]
    fold   = [int(a[1]) for a in amplitudes]
    nels   = sum(fold)

    diff = set(labels)-expected
    if diff:
        raise ValueError('Unexpected amplitude label(s) (%s), expected %s.' % (diff, expected))

    # all combinations
    all_permutations = product(permutations(range(1, nels+1)), permutations(range(1, nels+1)))
    components = [configuration_split(occ, vac, fold, labels) for occ, vac in all_permutations]

    # filter repeating components
    result = set(components)

    if verbose:
        print('%4d component(s) of kind %s' % (len(result), kind), file=sys.__stdout__)

    return result

####################################################################################################

def standardize_component(component):

    new = sorted(component, key=sort_order)

    occ = [orb for mult in new for orb in mult[1]]
    vac = [orb for mult in new for orb in mult[2]]

    ln = len(occ)

    # resort
    permuts = 0
    for orbs in (occ, vac):
        while True:
            current = 0
            for k in range(ln-1):
                if orbs[k] > orbs[k+1]:
                    current += 1
                    orbs[k], orbs[k+1] = orbs[k+1], orbs[k]

            # sorted, no permutations carried on current iteration
            if not current:
                break
            permuts += current

    sign = '+' if (-1)**permuts > 0 else '-'

    return sign, new

####################################################################################################

def generate_component_codeline(raw, scheme='alpha'):
    sign, component = standardize_component(raw)

    if scheme == 'alpha':
        # use either i,j,k,l,m,n,o,p / a,b,c,d,e,f,g,h
        omap = dict(zip(range(1,9), 'ijklmnop'))
        vmap = dict(zip(range(1,9), 'abcdefgh'))
    elif scheme == 'alphanumerical':
        # or i1,...,i8 / a1,...,a8
        omap = dict(zip(range(1,9), ['i%s' % k for k in range(1, 9)]))
        vmap = dict(zip(range(1,9), ['a%s' % k for k in range(1, 9)]))
    else:
        msg = 'Expected <scheme> to be one of (alpha, alphanumerical), while got %s.' % scheme
        raise ValueError(msg)

    row = []
    for mult in component:
        label, occ, vac = mult
        row.append('%s%s(%s,%s)' % (label, len(occ),
                                    ','.join([omap[o] for o in occ]),
                                    ','.join([vmap[v] for v in vac])))
    return sign, '*'.join(row)

####################################################################################################

def prepare_components(coefficient, verbose=False):
    components = []
    for component in generate_components(coefficient):
        components.append(generate_component_codeline(component, 'alphanumerical'))

    components = sorted(components, key=lambda x: x[1])
    if verbose:
        for sign, line in components:
            print(sign+line, file=sys.__stdout__)

    return components

####################################################################################################

def generate_coefficient_function(coefficient, prefix='cc_', indent='    '):

    def csort_order(mult):
        o1, o2 = ord(mult[0])-ord('a'), int(mult[1])
        return 1000*o1-10*o2

    rtype = 'real(kind=rglu)'
    itype = 'integer(kind=iglu)'

    components = prepare_components(coefficient, False)
    order      = sum([int(a[1]) for a in coefficient.split('*')])
    hamp       = ''.join(sorted(coefficient.split('*'), key=csort_order))
    name       = '%sc%s_%s' % (prefix, int(order), hamp)
    occupied   = ['i%s' % k for k in range(1, order+1)]
    vacant     = ['a%s' % k for k in range(1, order+1)]
    indices    = ','.join(occupied + vacant)

    print('%s%s %s %s(%s) result(ret)' % (indent, rtype, 'function', name, 'arr'))
    print('%s%s' %     (indent, 'implicit none'))
    print()
    print('%s%s, %s' % (indent, itype, 'intent(in) :: arr(:)'))
    print('%s%s  %s :: %s' %     (indent, itype, ' '*len('intent(in)'), indices))
    print()
    print()
    print('%s%s' % (indent, 'ret=0'))
    print()

    amps = coefficient.split('*')

    if 't1' in amps:
        print(indent + 'if (.NOT.pattern(1)) return')

    if 't2' in amps:
        print(indent + 'if (.NOT.pattern(2)) return')

    if 't3' in amps:
        print(indent + 'if (.NOT.pattern(3)) return')

    if 't1' in amps or 't2' in amps or 't3' in amps:
        print()

    print(indent + '; '.join(['i%s=arr(%s)' % (k, k) for k in range(1, order+1)]))
    print(indent + '; '.join(['a%s=arr(%s)' % (k, order+k) for k in range(1, order+1)]))
    print()
    if len(components) > 1:
        print('%s%s%s%s&' % (indent, 'ret=ret', components[0][0], components[0][1]))
        for sign, line in components[1:-1]:
            print('%s%s%s%s&' % (indent, ' '*len('ret=ret'), sign, line))
        print('%s%s%s%s' % (indent, ' '*len('ret=ret'), components[-1][0], components[-1][1]))
    else:
        print('%s%s%s%s' % (indent, 'ret=ret', components[0][0], components[0][1]))
    print()
    print('%s%s' % (indent, 'return'))
    print('%send function %s' % (indent, name))
    print()

####################################################################################################

def generate_include(pairs):

    print()
    for k, (prefix, coefficient) in enumerate(pairs):
        if k != 0:
            print('!   '+'~'*72+'   !\n')
        generate_coefficient_function(coefficient, prefix)

####################################################################################################

def dump_cc_configurations(file):

    svstdout, sys.stdout = sys.stdout, open(file, 'w')
    pairs = [
             ('cc_', 't1'),
             ('cc_', 't2'), ('cc_', 't1*t1'),
             ('cc_', 't3'), ('cc_', 't2*t1'), ('cc_', 't1*t1*t1'),
             ('cc_', 't3*t1'), ('cc_', 't2*t2'), ('cc_', 't2*t1*t1'), ('cc_', 't1*t1*t1*t1'),
            ]
    generate_include(pairs)
    sys.stdout = svstdout

####################################################################################################

def dump_lr_configurations(file):

    svstdout, sys.stdout = sys.stdout, open(file, 'w')
    pairs = [
             ('lr_', 'r1'),
             ('lr_', 'r2'), ('lr_', 'r1*t1'),
             ('lr_', 'r2*t1'), ('lr_', 'r1*t2'), ('lr_', 'r1*t1*t1'),
             ('lr_', 'r2*t2'), ('lr_', 'r2*t1*t1'), ('lr_', 'r1*t2*t1'), ('lr_', 'r1*t1*t1*t1'),
            ]
    generate_include(pairs)
    sys.stdout = svstdout

####################################################################################################

if __name__ == '__main__':

    dump_cc_configurations('cc_configurations.inc')
    dump_lr_configurations('lr_configurations.inc')
