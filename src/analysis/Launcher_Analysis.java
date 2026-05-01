package analysis;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.HashMap;

import optimisation.MenB_RMP_NM_Optimistion;
import sim.Simulation_ClusterModelTransmission;

public class Launcher_Analysis {

	public static void main(String[] args) throws IOException {
		if ("-opt".equals(args[0]) ? args.length < 3 : args.length < 2) {
			System.out.println(
					"Usage: java -jar Analyse_MenB_Vaccine.jar BASEDIR_SIM PATH_REGION_MAPPING PATH_GRP_SIZE <-flag=BINARY_FLAG_OPTIONS> <-printProgress=TF>"
							+ "\n  or java -jar Analyse_MenB_Vaccine.jar -opt BASEDIR_SIM SEED_DIR_NAME <-printProgress=TF> <-optType=OPT_TYPE>");
			System.exit(0);
		} else {
			if ("-opt".equals(args[0])) {
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

			} else {
				// Time trend
				Analysis_PostSim_ExtractTimeTrends analyse = new Analysis_PostSim_ExtractTimeTrends();
				analyse.analyse(args);

				// Infection history for PID
				Analysis_PostSim_ExtractInfectionHistory analyseInfHist = new Analysis_PostSim_ExtractInfectionHistory(
						args);

				// PID
				int[] incl_start_grps = new int[] { 5, 6, 7, 8, 9 }; // Indigenous female
				int[] sample_time = new int[] { 7300, 7665, 8030, 8395, 8760, 9125, 9490, 9855, 10220, 10585, 10950 };

				int max_exposure = 120; // Assume won't develop PID after 4 months
				double[] event_prob_by_inf_count = new double[] { 0.14, 0.17 };
				int[] inf_count_range = new int[] { 0, 1 };

				HashMap<String, double[]> resmap = analyseInfHist.analyse(incl_start_grps, sample_time, max_exposure,
						event_prob_by_inf_count, inf_count_range);

				if (!resmap.isEmpty()) {
					File baseDir = new File(args[0]);
					PrintWriter pWri = new PrintWriter(new File(baseDir, "PID_Event_Count.csv"));
					String[] zipEntNames = resmap.keySet().toArray(new String[0]);
					Arrays.sort(zipEntNames);

					pWri.print("Time");
					for (int t = 0; t < sample_time.length; t++) {
						pWri.print(',');
						pWri.print(sample_time[t]);
					}
					pWri.println();

					for (String zName : zipEntNames) {
						pWri.print(zName.replace(',', '_'));
						try {
							double[] ent = resmap.get(zName);
							for (int t = 0; t < ent.length; t++) {
								pWri.print(',');
								pWri.print(ent[t]);
							}
						} catch (NullPointerException ex) {
							ex.printStackTrace(System.err);
						}
						pWri.println();
					}

					pWri.close();

				}

			}

		}

	}

}
