"""Narrow vendor contract and unchanged background delegation tests."""
import unittest
from unittest.mock import patch
import numpy as np
from ape_c_lpsd.analysis import vendor_statistics as adapter


class DummyFit:
    @staticmethod
    def process_iteration(args):
        _, x, y, config = args
        assert len(x) == len(y) and config["buffer"] == 50
        return 0, np.array([x[0], x[-1], y[0], y[-1]]), [0., .5, 2.], 1.5

    @staticmethod
    def bkg_model(x, knots):
        return np.full_like(x, -4.)


class DummyStats:
    @staticmethod
    def log_likelihood(mu, Y, bkg, peak_norm, peak_shape, model_args):
        return float(-np.sum((Y-np.log(np.exp(bkg)+mu*peak_norm*peak_shape))**2))


class IntegrationTests(unittest.TestCase):
    def test_background_preserved(self):
        f = np.linspace(100, 101, 1100)
        p = np.exp(np.full(1100, -4.))
        frozen = p.copy()
        with patch.object(adapter, "load_external", return_value=(DummyFit,DummyStats,{})):
            result = adapter.fit_window(f, p, "/unused")
        np.testing.assert_array_equal(p, frozen)
        self.assertTrue(result["vendor_quality_pass"])

    def test_multi_observation_frequency_dependent_norm(self):
        Y = np.log(np.full((2, 4), 2.))
        bkg = np.log(np.ones((2, 4)))
        norm = np.array([[.1,.2,.3,.4],[.5,.4,.3,.2]])
        shape = np.array([.1,.2,.3,.4])
        args = np.array([[0.,1.,0.],[0.,1.,0.]])
        with patch.object(adapter, "load_external", return_value=(DummyFit,DummyStats,{})):
            observed = adapter.log_likelihood(.3,Y,bkg,norm,shape,args,"/unused")
        self.assertAlmostEqual(observed,DummyStats.log_likelihood(.3,Y,bkg,norm,shape,args))
        with self.assertRaises(ValueError):
            adapter.log_likelihood(.3,Y,bkg,norm,np.ones(4),args,"/unused")

    def test_signal_template_factorization(self):
        power=np.array([[.1,.3,0.],[.4,.2,0.]])
        norm,shape=adapter.factor_power_template(power)
        np.testing.assert_allclose(norm*shape[None,:],power,atol=1e-15)
        self.assertAlmostEqual(shape.sum(),1.)


if __name__ == "__main__":
    unittest.main()
