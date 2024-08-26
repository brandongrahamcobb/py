from datetime import datetime, timedelta

class CooldownManager:
    def __init__(self, cooldown_seconds=30):
        self.cooldowns = {}
        self.cooldown_time = timedelta(seconds=cooldown_seconds)

    def check_cooldown(self, user_id):
        now = datetime.now()
        if user_id in self.cooldowns:
            last_message_time = self.cooldowns[user_id]
            if now - last_message_time < self.cooldown_time:
                return False, self.cooldown_time - (now - last_message_time)
        self.cooldowns[user_id] = now
        return True, None
