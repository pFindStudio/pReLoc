__version__ = '0.26.0+cpu'
git_version = '336d36e8db990a905498c73933e35231876e28bc'
from torchvision.extension import _check_cuda_version
if _check_cuda_version() > 0:
    cuda = _check_cuda_version()
