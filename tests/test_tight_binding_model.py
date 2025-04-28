import pytest
import numpy as np
from pyamtb.tight_binding_model import TightBindingModel
from pyamtb.parameters import Parameters

def test_model_initialization():
    """Test TightBindingModel initialization"""
    params = Parameters()
    model = TightBindingModel(params)
    assert isinstance(model, TightBindingModel)
    assert hasattr(model, 'parameters')
    assert hasattr(model, 'hamiltonian')

def test_hamiltonian_shape():
    """Test Hamiltonian matrix shape"""
    params = Parameters()
    model = TightBindingModel(params)
    k = np.array([0.0, 0.0, 0.0])
    H = model.hamiltonian(k)
    assert isinstance(H, np.ndarray)
    assert H.shape[0] == H.shape[1]  # Hamiltonian should be square matrix

def test_energy_bands():
    """Test energy bands calculation"""
    params = Parameters()
    model = TightBindingModel(params)
    k_path = np.array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]])
    bands = model.calculate_bands(k_path)
    assert isinstance(bands, np.ndarray)
    assert len(bands.shape) == 2  # Should be 2D array (k-points x bands) 