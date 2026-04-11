package analysis;

import java.io.IOException;

import optimisation.MenB_RMP_NM_Optimistion;

public class Launcher_Analysis {

	public static void main(String[] args) throws IOException {
		if (args.length < 3) {
			System.out.println(
					"Usage: java -jar Analyse_MenB_Vaccine.jar BASEDIR_SIM PATH_REGION_MAPPING PATH_GRP_SIZE <-flag=BINARY_FLAG_OPTIONS> <-suppressOutput=TF>"
							+ "\n  or java -jar Analyse_MenB_Vaccine.jar -opt BASEDIR_SIM SEED_DIR_NAME <-suppressOutput=TF>");
			System.exit(0);
		} else {
			if ("-opt".equals(args[0])) {										
				MenB_RMP_NM_Optimistion opt = new MenB_RMP_NM_Optimistion(args[1], args[2]);				
				for(int i = 2; i < args.length; i++) {
					if(args[i].startsWith("-suppressOutput")) {											
						opt.setPrintProgress(Boolean.parseBoolean(args[i].split("=")[1]));
					}
				}
				opt.runOptimisation();				

			} else {
				Analysis_PostSim_ExtractResults analyse = new Analysis_PostSim_ExtractResults();
				analyse.analyse(args);
			}

		}

	}

}
