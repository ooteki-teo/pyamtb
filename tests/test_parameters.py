import pytest
import numpy as np
from pyamtb.parameters import Parameters

def test_parameters_initialization():
    """Test Parameters class initialization"""
    params = Parameters()
    assert isinstance(params, Parameters)
    assert hasattr(params, 'lattice_constant')
    assert hasattr(params, 'hopping_parameters')

def test_parameters_from_toml():
    """Test loading parameters from TOML file"""
    # Create a temporary TOML file
    import tempfile
    import tomlkit
    
    test_params = {
        'lattice_constant': 1.0,
        'hopping_parameters': {
            't1': 1.0,
            't2': 0.5
        }
    }
    
    with tempfile.NamedTemporaryFile(suffix='.toml', mode='w', delete=False) as f:
        tomlkit.dump(test_params, f)
        f.flush()
        
        params = Parameters.from_toml(f.name)
        assert params.lattice_constant == test_params['lattice_constant']
        assert params.hopping_parameters['t1'] == test_params['hopping_parameters']['t1']
        assert params.hopping_parameters['t2'] == test_params['hopping_parameters']['t2'] 