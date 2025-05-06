import numpy as np

class StateTracker():
    def __init__(self):
        # Stores the current state, increments when state is updated
        self.currentStateCounter = 0
        self.id = np.random.random()
        self.currentState = {}

    def updateState(self, newState):
        # Reject if updates are from different tracker
        if newState['id'] != self.id:
            return False
        
        # Track if any modifications were made
        modified = False

        # Extract counter
        counter = newState.pop('counter')

        # We always take previously undefined keys
        for k, v in newState.items():
            if k not in self.currentState:
                self.currentState[k] = v
                modified = True

        # We only accept conflicts for new states
        if counter >= self.currentStateCounter:
            conflicts = {
                k: newState[k] for k in newState.keys() 
                if self.currentState[k] != newState[k]
            }

            if len(conflicts) != 0:
                modified = True

            for k, v in conflicts.items():
                self.currentState[k] = v
        
        if modified:
            self.currentStateCounter += 1

        
        if modified:
            return True
        else:
            return False
        
    def getState(self):
        # Never return the original object, deepcopy shouldnt be 
        # neccessary as dict is not nested
        state = self.currentState.copy()
        state['counter'] = self.currentStateCounter
        state['id'] = self.id
        return state
