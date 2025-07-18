#print("fingal imported")
from .tools import *
from .datahandling import *
from .tempinversion import *
# #from .fieldERT import *
from .inversionsIP import *
from .inversionsIP2 import *
from .inversionsERT import *
from .inversionsERTWithRobinCondition import ERTInversionH1WithRobinCondition
from .datamapping import mapToDomain
from .synthetics import IPSynthetic
from .sensitivity import ERTSensitivity
from .geometry import MeshWithTopgraphy
