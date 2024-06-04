CHANNEL_CONF = {
    'CE89324': {
        'Name': 'Rio Dell - Hwy 101/Painter St. Overpass',
        'Transverse': {
            'inputs':  [3,17,20],
            'outputs': [7,9,4]
        },
        'Longitudinal': {
            'inputs':  [15,1,18],
            'outputs': [11]
        }
    },
    'CE01336': {
        'Name': 'El Centro - Hwy8/Meloland Overpass',
        'Transverse 1': {
            'inputs':  [2],
            'outputs': [5,7,9]
        },
        'Transverse 2': {
            'inputs':  [11,2,26],
            'outputs': [5,7,9]
        },
        'Longitudinal': {
            'inputs':  [12,4,25],
            'outputs': [27,8]
        },
    },
    'CE54730': {
        'Name': 'Lake Crowley - Hwy 395 Bridge',
        'Transverse 1': {
            'inputs':  [4],
            'outputs': [6,7,9]
        },
        'Transverse 2': {
            'inputs':  [6,4,9],
            'outputs': [7]
        },
        'Transverse 3': {
            'inputs':  [4],
            'outputs': [7]
        },
        'Longitudinal': {
            'inputs':  [5],
            'outputs': [8]
        },
    },
    'CE33742': {
        'Name': 'Ridgecrest - Hwy 395/Brown Road Bridge',
        'Transverse': {
            'inputs':  [4],
            'outputs': [6,7,9]
        },
    },
    'CE13795': {
        'Name': 'Capistrano Beach - I5/Via Calif. Bridge',
        'Transverse': {
            'inputs':  [4],
            'outputs': [10,7]
        },
    },
    'CE58658': {
        'Name': 'Hayward - Hwy 580/238 Interchange Bridge',
        'Transverse': {
            'inputs':  [25,2,7,18],
            'outputs': [23,13,15,20]
        },
        'Longitudinal': {
            'inputs':  [3,6,17],
            'outputs': [12,14,19]  # For SISO, try 6,14 in addition to default 3,12
        },
    },
    'CE23631': {
        'Name': 'San Bernardino - I10/215 Interchange Br',
        'Transverse 1': {
            'inputs':  [6],
            'outputs': [7,8],
            'description': 'Bent 3'
        },
        'Transverse 2': {
            'inputs':  [24],
            'outputs': [19,20],
            'description': 'Bent 8'
        },
        'Longitudinal 1': {
            'inputs':  [4],
            'outputs': [10],
            'description': 'Bent 3'
        },
        'Longitudinal 2': {
            'inputs':  [22],
            'outputs': [17,18],
            'description': 'Bent 8'
        },
    },
    'CE14406': {
        'Name': "Los Angeles - Vincent Thomas Bridge",
        'Transverse 1': {
            'inputs':  [1,9,24],
            'outputs': [2,5,7]
        },
        'Transverse 2': {
            'inputs':  [1,9,24],
            'outputs': [2,4,5,6,7],
            'description': 'Dense'
        },
        'Vertical': {
            'inputs':  [14,19,26],
            'outputs': [16,18,22]
        },
    },
}