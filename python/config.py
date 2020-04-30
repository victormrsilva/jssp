class Config:
    def __init__(self, file):
        self.conf = {}
        with open(file) as fp:
            for line in fp.read().splitlines():
                if line.startswith('#') or line.startswith('[') or len(line) < 1:
                    continue
                key, val = line.strip().split('=')
                key = key.strip().lower()
                val = val.strip()
                if val.isdigit():  # if is integer
                    self.conf[key] = int(val)
                elif val.replace('.', '', 1).isdigit():  # if is decimal
                    self.conf[key] = float(val)
                else:  # if is integer
                    self.conf[key] = val

    def get_property(self, prop):
        return self.conf.get(prop, None)

