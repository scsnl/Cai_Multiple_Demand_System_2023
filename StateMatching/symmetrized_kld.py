import numpy as np
from numpy.linalg import slogdet, inv


class SymmetrizedKullbackLeiblerDivergenceGaussians(object):
    name = "sKLD_Gaussians"

    """ 
    sKLD = 0.5 ( KL(p || q) + KL(q || p) )
    """
    def __init__(self, cov1, cov2, mean1, mean2):
        """

        :param cov1: Covariance of Gaussian p1()
        :param cov2: Covariance of Gaussian p2()
        :param mean1: Mean of Gaussian p1()
        :param mean2: Mean of Gaussian p2()
        """
        if cov1.shape[0] != cov2.shape[0]:
            raise ValueError("dim mismatch")
        self.dim = cov1.shape[0]

        is_pd(cov2)
        is_pd(cov1)

        self.cov1 = cov1
        self.cov2 = cov2
        self.mean1 = mean1
        self.mean2 = mean2

    def evaluate(self):
        kld12 = self._kld_pq(cov_p=self.cov1, cov_q=self.cov2, mean_p=self.mean1, mean_q=self.mean2)
        kld21 = self._kld_pq(cov_p=self.cov2, cov_q=self.cov1, mean_p=self.mean2, mean_q=self.mean1)
        return 0.5 * (kld12 + kld21)

    def _kld_pq(self, cov_p, cov_q, mean_p, mean_q):
        inv_cov_p = inv(cov_p)
        term1 = slogdet(cov_p)[1] - slogdet(cov_q)[1]
        term2 = np.trace(inv_cov_p @ cov_q)
        term3 = (mean_p - mean_q) @ inv_cov_p @ (mean_p - mean_q)
        return 0.5 * (-self.dim + term1 + term2 + term3)


def is_pd(x):
    try:
        np.linalg.cholesky(x)
        return True
    except np.linalg.linalg.LinAlgError:
        return False


