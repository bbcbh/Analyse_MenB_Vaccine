package analysis;

import java.io.File;
import java.io.FileFilter;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.Properties;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import sim.Runnable_MetaPopulation_MultiTransmission;
import sim.SimulationInterface;
import util.Util_7Z_CSV_Entry_Extract_Callable;

public class Analysis_PostSim_ExtractInfectionHistory {

	private Pattern res_dir_format = Pattern.compile("Seed_List_(\\d+)");
	private Pattern infect_hist_7z_format = Pattern.compile("InfectHist_(.*)\\.csv\\.7z");
	private Pattern infect_hist_key = Pattern.compile("\\[(.*),(\\d+)\\]InfectHist_(-?\\d+)_(-?\\d+).csv");

	// private final int INFECTION_HIST_CLEAR_NATURAL_RECOVERY =
	// Runnable_MetaPopulation_MultiTransmission.INFECTION_HIST_CLEAR_NATURAL_RECOVERY;
	// private final int INFECTION_HIST_CLEAR_TREATMENT =
	// Runnable_MetaPopulation_MultiTransmission.INFECTION_HIST_CLEAR_TREATMENT;
	// From Runnable_MetaPopulation_MultiTransmission

	private File basedir;

	private Properties prop;

	private HashMap<Long, HashMap<Integer, int[]>> map_indiv_stat;
	private HashMap<String, ArrayList<String[]>> map_infhist_lines;

	public Analysis_PostSim_ExtractInfectionHistory(String[] args) {
		basedir = new File(args[0]);

		// K = CMAP_SEED, V = ID, stat
		// From Runnable_MetaPopulation_MultiTransmission
		map_indiv_stat = new HashMap<>();
		map_infhist_lines = new HashMap<>();

	}

	public HashMap<String, double[]> analyse() throws IOException {

		int[] incl_start_grps = new int[] { 5, 6, 7, 8, 9 }; // Indigenous female
		int[] sample_time = new int[] { 7300, 7665, 8030, 8395, 8760, 9125, 9490, 9855, 10220, 10585, 10950 };

		int max_exposure = 120; // Assume won't develop PID after 4 months
		double[] event_prob_by_inf_count = new double[] { 0.14, 0.17 };
		int[] inf_count_range = new int[] { 0, 1 };

		// Code start
		if (map_infhist_lines.isEmpty()) {
			loadInfHistMap();
		}
		HashMap<String, double[]> map_event_prob_by_file = new HashMap<>();

		for (Entry<String, ArrayList<String[]>> ent : map_infhist_lines.entrySet()) {
			Matcher m_map_infhist = infect_hist_key.matcher(ent.getKey());
			if (m_map_infhist.matches()) {
				Long cMap = Long.valueOf(m_map_infhist.group(3));
				HashMap<Integer, int[]> indivMap = map_indiv_stat.get(cMap);
				double[] event_prob_sum = new double[sample_time.length];
				map_event_prob_by_file.put(ent.getKey(), event_prob_sum);

				ArrayList<String[]> inf_hist_rows = ent.getValue();

				for (int r = 1; r < inf_hist_rows.size(); r++) {
					String[] row = inf_hist_rows.get(r);
					Integer id = Integer.parseInt(row[0]);
					int[] indiv_ent = indivMap.get(id);
					int start_grp = indiv_ent[Runnable_MetaPopulation_MultiTransmission.INDIV_MAP_ENTER_GRP];
					if (Arrays.binarySearch(incl_start_grps, start_grp) >= 0) {
						int sample_end_point = 0; // i.e. to be include in the sample entry
						ArrayList<Double> cumul_prob = new ArrayList<>(); // P(0), P(1) etc... within current window
						cumul_prob.add(Double.valueOf(1)); // P(0) = 1
						int inf_count = 0;
						for (int inf_pt = 2; inf_pt < row.length; inf_pt += 3) { //
							int exposure_start = Integer.parseInt(row[inf_pt]);

							// Stop checking if inf start after sample period
							if (exposure_start >= sample_time[sample_time.length - 1]) {
								break;
							}

							int exposure_end = exposure_start + max_exposure;
							// For case where individual leave population before recovery
							if (inf_pt + 1 >= row.length) {
								exposure_end = Math.min(exposure_end,
										indiv_ent[Runnable_MetaPopulation_MultiTransmission.INDIV_MAP_EXIT_POP_AT]);
							} else {
								exposure_end = Math.min(exposure_end, Integer.parseInt(row[inf_pt + 1]));
							}

							if (exposure_end >= sample_time[0]) {								
								exposure_start = Math.max(exposure_start, sample_time[0]);								
								int event_count_pt = Arrays.binarySearch(inf_count_range, inf_count);
								if (event_count_pt < 0) {
									event_count_pt = ~event_count_pt;
								}
								double current_event_prob = event_prob_by_inf_count[Math
										.min(event_prob_by_inf_count.length - 1, event_count_pt)];

								while (exposure_start > sample_time[sample_end_point]) {
									sample_end_point = updateEventProbSum(event_prob_sum, sample_end_point, cumul_prob);
								}

								// Check for reduced exposure
								if (exposure_end < exposure_start + max_exposure) {
									double pDayExp = 1 - Math.exp(Math.log(1 - current_event_prob) / max_exposure);
									current_event_prob = 1 - Math.pow(1 - pDayExp, exposure_end - exposure_start);
								}

								// Update cumul_prod
								cumul_prob.add(cumul_prob.get(cumul_prob.size() - 1) * current_event_prob);
								for (int pI = cumul_prob.size() - 2; pI >= 1; pI--) {
									cumul_prob.set(pI, cumul_prob.get(pI) * (1 - current_event_prob)
											+ cumul_prob.get(pI - 1) * current_event_prob);
								}
								cumul_prob.set(0, cumul_prob.get(0) * (1 - current_event_prob));

								int expoure_to_next_window = Math.max(exposure_end - sample_time[sample_end_point], 0);

								if (expoure_to_next_window > 0) {
									sample_end_point = updateEventProbSum(event_prob_sum, sample_end_point, cumul_prob);
									double pDayExp = 1 - Math.exp(Math.log(1 - current_event_prob) / max_exposure);
									current_event_prob = 1 - Math.pow(1 - pDayExp, expoure_to_next_window);
									cumul_prob.add(current_event_prob);
									cumul_prob.set(0, cumul_prob.get(0) * (1 - current_event_prob));
								}
							}
							inf_count++;
						}

						sample_end_point = updateEventProbSum(event_prob_sum, sample_end_point, cumul_prob);

					} // End of if (Arrays.binarySearch(incl_start_grps, start_grp) >= 0) {

				} // End of for (int r = 1; r < inf_hist_rows.size(); r++) {
			}
		}
		return map_event_prob_by_file;

	}

	private int updateEventProbSum(double[] event_prob_sum, int sample_end_point, ArrayList<Double> cumul_prob) {
		// Update event_prob_sum
		double weighted_mean = 0;
		for (int i = 1; i < cumul_prob.size(); i++) {
			weighted_mean += i * cumul_prob.get(i);
		}
		event_prob_sum[sample_end_point] += weighted_mean;
		// Move event window
		sample_end_point++;
		// Clear stat
		cumul_prob.clear();
		cumul_prob.add(Double.valueOf(1)); // P(0) = 1
		return sample_end_point;
	}

	private void loadInfHistMap() throws IOException {

		// Reading of PROP file
		File propFile = new File(basedir, SimulationInterface.FILENAME_PROP);
		FileInputStream fIS = new FileInputStream(propFile);
		prop = new Properties();
		prop.loadFromXML(fIS);
		fIS.close();

		File[] res_dirs = basedir.listFiles(new FileFilter() {
			@Override
			public boolean accept(File pathname) {
				return pathname.isDirectory() && res_dir_format.matcher(pathname.getName()).matches();
			}
		});

		for (File res_dir : res_dirs) {
			File[] file_infhist_7z = res_dir.listFiles(new FileFilter() {
				@Override
				public boolean accept(File pathname) {
					return infect_hist_7z_format.matcher(pathname.getName()).matches();
				}
			});

			for (File zip : file_infhist_7z) {
				map_infhist_lines = Util_7Z_CSV_Entry_Extract_Callable.extractedLinesFrom7Zip(zip, map_infhist_lines);
			}
		}

		for (Entry<String, ArrayList<String[]>> ent : map_infhist_lines.entrySet()) {
			Matcher m_map_infhist = infect_hist_key.matcher(ent.getKey());
			if (m_map_infhist.matches()) {
				Long cMap = Long.valueOf(m_map_infhist.group(3));

				HashMap<Integer, int[]> indivMap = map_indiv_stat.get(cMap);
				if (indivMap == null) {
					indivMap = new HashMap<>();
					map_indiv_stat.put(cMap, indivMap);

					File demo_dir = new File(prop.getProperty("PROP_CONTACT_MAP_LOC"));
					demo_dir = new File(demo_dir, String.format("Demographic_%d", cMap));
					String[] pop_stat_line = Util_7Z_CSV_Entry_Extract_Callable
							.extracted_lines_from_text(new File(demo_dir, String.format("POP_STAT_%d.csv", cMap)));
					// ID,GRP,ENTER_POP_AGE,ENTER_POP_AT,EXIT_POP_AT,HOME_LOC
					for (int i = 1; i < pop_stat_line.length; i++) {
						String[] lineEnt = pop_stat_line[i].split(",");
						int[] indiv_ent = new int[Runnable_MetaPopulation_MultiTransmission.LENGTH_INDIV_MAP];
						indivMap.put(Integer.valueOf(lineEnt[0]), indiv_ent);
						indiv_ent[Runnable_MetaPopulation_MultiTransmission.INDIV_MAP_ENTER_GRP] = Integer
								.parseInt(lineEnt[1]);
						indiv_ent[Runnable_MetaPopulation_MultiTransmission.INDIV_MAP_ENTER_POP_AGE] = Integer
								.parseInt(lineEnt[2]);
						indiv_ent[Runnable_MetaPopulation_MultiTransmission.INDIV_MAP_ENTER_POP_AT] = Integer
								.parseInt(lineEnt[3]);
						indiv_ent[Runnable_MetaPopulation_MultiTransmission.INDIV_MAP_EXIT_POP_AT] = Integer
								.parseInt(lineEnt[4]);
						indiv_ent[Runnable_MetaPopulation_MultiTransmission.INDIV_MAP_HOME_LOC] = Integer
								.parseInt(lineEnt[5]);
					}
				}

			} else {
				System.err.printf("Warning! Illformed zip file entry %s. Entry ignored.\n", ent.getKey());
			}
		}
	}

}
