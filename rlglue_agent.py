import sys
sys.path.append('/Users/tomii/Source/Python/dqn_chainer/dqn_chainer')

from rlglue.agent import AgentLoader as AgentLoader
from dqn_agent import dqn_agent
from optparse import OptionParser

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option(
        '-t','--test', 
        dest='test', action='store_true',default=False,
        help="test mode")
    parser.add_option(
        '-p','--param_file', 
        dest='param_file', action='store', type='string', default='./rlglue_param.json',
        help="json file of simulation parameters")
    parser.add_option(
        '-d','--dst', 
        dest='savepath', action='store', type='string', default='./result/data',
        help="Save data path.")
    (options, args) = parser.parse_args()
    print 'options', options
  
    if not options.test:
        AgentLoader.loadAgent(dqn_agent())
    else:
        from rlglue.types import Observation
        
        print 'dqn_agent test'
        
        imgSize = 64
        obs = Observation(numDoubles = imgSize ** 2) 
        objAgent = dqn_agent()
        objAgent.agent_init(str(imgSize))
        objAgent.agent_start(obs)
        while True:
          objAgent.agent_step(0.0, obs)
