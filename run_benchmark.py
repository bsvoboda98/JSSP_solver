import slurminade
import subprocess

slurminade.update_default_configuration(
    partition="alg",
    constraint="alggen05",
    mail_user="b.svoboda@tu-braunschweig.de",
    mail_type="FAIL",
    exclusive=True
)

slurminade.set_dispatch_limit(1_000)
@slurminade.slurmify()
def run(arguments):
    #needs to be adapted later on
    executable_path = '/ibr/home/svoboda/bachelorarbeit-bennet-svoboda/JSSP_solver/cmake-build-release/JSSP_solver'

    command = [executable_path] + arguments
    try:
        result = subprocess.run(command, check=True, capture_output=False, cwd='/ibr/home/svoboda/bachelorarbeit-bennet-svoboda/JSSP_solver/cmake-build-release', stderr=subprocess.PIPE)

    except subprocess.CalledProcessError as e:
        print("Error:", e)
        print("Return code:", e.returncode)
        print("Output:", e.output.decode())
        print("Error output:", e.stderr.decode())




if __name__ == "__main__":




    instanceList = ["dmu37-40", "dmu41-45", "dmu46-50", "dmu51-55", "dmu56-60", "dmu61-65", "dmu66-70", "la01-05", "la06-10", "ta31-35", "ta36-40", "ta41-45"]
    for instance in instanceList:

        ga_arguments = ["-ga", "-g", "20", "0", "0.7", "0.7", "0.5", "-i", instance, "-f", "cc", "ori", "rs", "-t", "1200", "-bm"]
        run.distribute(ga_arguments)

        ga_arguments = ["-ga", "-g", "20", "0", "0.7", "0.7", "0.4", "-i", instance, "-f", "cc", "ori", "rs", "-t", "1200", "-bm"]
        run.distribute(ga_arguments)

        or_arguments = ["-or", "-i", instance, "-t", "1200", "-bm"]
