import logging

import numpy as np

logger = logging.getLogger(__name__)


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

        # We always take previously undefined keys. A None ("nothing selected") must
        # not claim a key: the frontend sends null for every unset field, so the first
        # cell to report would otherwise own every key, and a later first-time value
        # carrying the same counter (e.g. the Scan Table's default row) would be
        # treated as a stale conflict below and dropped.
        for k, v in newState.items():
            if k == 'id' or v is None:
                continue
            if k not in self.currentState:
                self.currentState[k] = v
                modified = True

        conflicts = {
            k: newState[k] for k in newState.keys()
            if k in self.currentState and self.currentState[k] != newState[k]
        }

        # We only accept conflicts for new states
        if counter >= self.currentStateCounter:
            if len(conflicts) != 0:
                modified = True

            for k, v in conflicts.items():
                self.currentState[k] = v
        elif conflicts:
            # A cell sent a change while it had not yet received the latest state
            # (another cell's update was in flight). The change is lost; log it so
            # "my click did nothing" reports can be checked against the server log.
            logger.debug(
                'StateTracker: ignored update with counter %s < %s for keys %s',
                counter, self.currentStateCounter, sorted(conflicts)
            )

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
