package analysis;

import java.io.IOException;

import optimisation.MenB_RMP_NM_Optimistion;
import sim.Simulation_ClusterModelTransmission;

public class Launcher_Analysis {

	public static void main(String[] args) throws IOException {
	
		String usageInfo = "Usage: java -jar Analyse_MenB_Vaccine.jar BASEDIR_SIM PATH_REGION_MAPPING PATH_GRP_SIZE <-flag=BINARY_FLAG_OPTIONS> <-printProgress=TF>"
				+ "\n  or java -jar Analyse_MenB_Vaccine.jar -opt BASEDIR_SIM SEED_DIR_NAME <-printProgress=TF> <-optType=OPT_TYPE>"
				+ "\n  or java -jar Analyse_MenB_Vaccine.jar -optBestFit BASEDIR_SIM CSV_FOLDER_NAME";
	
		boolean correctUsage = false;
		if (args.length > 0) {	
			if ("-opt".equals(args[0]) && args.length >= 3) {
				// Optimisation
				correctUsage = true;
				MenB_RMP_NM_Optimistion opt = new MenB_RMP_NM_Optimistion(args[1], args[2]);
				for (int i = 2; i < args.length; i++) {
					if (args[i].startsWith(Simulation_ClusterModelTransmission.LAUNCH_ARGS_PRINT_PROGRESS)) {
						opt.setObjFunc_PrintProgress(Boolean.parseBoolean(args[i].split("=")[1]));
					}
					if (args[i].startsWith("-optType")) {
						opt.setOptType(Integer.parseInt(args[i].split("=")[1]));
					}
				}
				opt.runOptimisation();
			} else if ("-optBestFit".equals(args[0]) && args.length >= 3) {
				// Optimisation result extraction
				correctUsage = true;
				Analysis_PostOptimsation.analysis(args[1], args[2]);														
			} else if (args.length >= 3) {
				// Time trend
				correctUsage = true;
				Analysis_PostSim_ExtractTimeTrends analyse = new Analysis_PostSim_ExtractTimeTrends();
				analyse.analyse(args);							
			}
		}
	
		if (!correctUsage) {
			System.out.println(usageInfo);
			System.exit(0);
		}
	
	}

}
