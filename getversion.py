
def define_env(env):

    @env.macro
    def spades_version():
        try:
            lines = open("VERSION").readlines()
            version = lines[0].strip()
            # FIXME: dirty hack for current VERSION file
            if version.find("dev") != -1:
                version = "3.15.5"
            return version
        except:
            return "3.15.5"

