'''XML utility functions'''


def convert_floats(dict_):
    '''convert to float where possible

    Note: original dict is modified in place
    '''
    for k, v in dict_.items():
        try:
            dict_[k] = float(v)
        except ValueError:
            pass
    return dict_


def read_placement_rules(tree):
    '''convert an placement_rules xml string to a dictionary of rules

    Format of dict is:
        {'L1_HAC':
            {'L1_HAC, axon, Layer_1':
                {'old_id': 'L1_HAC_axon_target',
                'segment_type': 'axon',
                'type': 'region_target',
                'y_max_fraction': '1.00',
                'y_max_layer': '1',
                'y_min_fraction': '0.00',
                'y_min_layer': '1'}},
    '''
    root = tree.getroot()
    assert root.tag == 'annotations'

    ret = {}

    for child in root.getchildren():
        if child.tag == 'placement':
            rule = child.get('rule')
            if 'axon' in rule:
                ret['axon'] = convert_floats(dict(child.items()))
            if 'dendrite' in rule:
                ret['dendrite'] = convert_floats(dict(child.items()))

    return ret


def update_rule(root, pattern, properties):
    '''Search from the root Element for the rule 'rule'
    and update all the attributes in the dictionary'''
    for child in root.getchildren():
        if child.tag == 'placement' and pattern in child.get('rule'):
            for key, value in properties.items():
                child.set(key, value)
            return
    raise RuntimeError('No rule found to update')
