import math

class Progress:

    def __init__(self, data, indicator):
        self.data = data
        if self.data:
            self.data['type'] = indicator
            if indicator == "percentage-indicator":
                self.data['percentage'] = 0
            self.data["status"] = "in-progress"
        return

    def update(self, percentage):
        if self.data:
            self.data["percentage"] = percentage
        return

    def update_n_of_m(self, m, n):
        if self.data:
            if m < n:
                total = n
            else:
                total = m
            self.data["percentage"] = math.trunc( (100 * m) / total )
        return

    def ready(self):
        if self.data:
            self.data["status"] = "ready"
        return

    def fail(self, err_msg):
        if self.data:
            self.data["status"] = "failed"
            self.data["error"] = err_msg
        return
