import unittest, numpy as np
from evaluate import landscape_metrics, finite, rank_correlation
class Metrics(unittest.TestCase):
 def test_s3(self):
  y=np.array([-1.,0.,1.]); p=np.arange(3)
  for x,mae,cmae in [(y+.5,.5,0),(y*.5,1/3,1/3),(y*0,2/3,2/3)]:
   v=landscape_metrics(y,x,p)
   self.assertAlmostEqual(v['ordinary_mae'],mae); self.assertAlmostEqual(v['cmae'],cmae)
 def test_position_weights(self):
  v=landscape_metrics(np.zeros(3),np.array([1.,3.,9.]),np.array(['a','a','b']))
  self.assertAlmostEqual(v['ordinary_mae'],5.5)
 def test_constant_and_finite(self):
  self.assertIsNone(rank_correlation([0,0,0],[1,2,3])); self.assertFalse(finite('nan')); self.assertFalse(finite('')); self.assertFalse(finite('inf'))
 def test_identities(self):
  v=landscape_metrics(np.array([-1.,0.,1.]),np.array([1.,0.,1.]),np.arange(3))
  self.assertAlmostEqual(v['total_mse'],v['mean_mse']+v['centered_mse'])
  self.assertAlmostEqual(v['centered_mse'],v['predicted_variance']+v['observed_variance']-2*v['covariance'])
if __name__=='__main__': unittest.main()
