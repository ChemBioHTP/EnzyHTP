import pytest
import lasso_extender as lasso
import numpy as np

def test_rotation_matrix_parallel():
    assert np.array_equiv(lasso.rotation_matrix_from_vectors(np.array([1, 1, 1]), np.array([2, 2, 2])), np.eye(3))

def test_rotation_matrix_aroundY():
    assert np.array_equal(lasso.rotation_matrix_from_vectors(np.array([0, 0, 2]), np.array([2, 0, 0])), np.array([[0, 0, 1], [0, 1, 0,], [-1, 0, 0]]))

def test_rotation_matrix_wrong_type():
    with pytest.raises(ValueError):
        assert lasso.rotation_matrix_from_vectors(0, 5)

def test_rotation_matrix_wrong_dim():
    with pytest.raises(ValueError):
        assert lasso.rotation_matrix_from_vectors(np.array([1, 1]), np.array([0, 0, 0, 4]))

print(type(np.array([0, 1, 1])))