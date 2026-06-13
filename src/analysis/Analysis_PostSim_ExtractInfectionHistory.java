package analysis;

import java.io.File;
import java.io.FileFilter;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.Properties;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.compress.archivers.sevenz.SevenZArchiveEntry;
import org.apache.commons.compress.archivers.sevenz.SevenZFile;

import sim.Runnable_MetaPopulation_MultiTransmission;
import sim.SimulationInterface;
import util.StaticMethods;

public class Analysis_PostSim_ExtractInfectionHistory {

	private Pattern res_dir_format = Pattern.compile("Seed_List.*_(\\d+)");
	private Pattern infect_hist_7z_format = Pattern.compile("InfectHist_(.*)\\.csv\\.7z");
	private Pattern infect_hist_key = Pattern.compile("\\[(.*),(\\d+)\\]InfectHist_(-?\\d+)_(-?\\d+).csv");

	private File basedir;

	private Properties prop;

	private HashMap<Long, HashMap<Integer, int[]>> map_indiv_stat;
	private HashMap<String, ArrayList<int[]>> map_infhist_lines;

	public static final Pattern patten_zipEnt = Pattern.compile("\\[Seed_List_(\\d+)\\.csv,(\\d+)\\].*");
	public static final Comparator<String> cmp_zipEnt = new Comparator<String>() {
		@Override
		public int compare(String o1, String o2) {
			int res = o1.compareTo(o2);
			Matcher m1 = patten_zipEnt.matcher(o1);
			Matcher m2 = patten_zipEnt.matcher(o2);

			if (m1.matches() && m2.matches()) {
				res = Integer.valueOf(m1.group(1)).compareTo(Integer.valueOf(m2.group(1)));
				if (res == 0) {
					res = Integer.valueOf(m1.group(2)).compareTo(Integer.valueOf(m2.group(2)));
				}

			}

			return res;
		}
	};

	public Analysis_PostSim_ExtractInfectionHistory(String[] args) {
		basedir = new File(args[0]);

		// K = CMAP_SEED, V = ID, stat
		// From Runnable_MetaPopulation_MultiTransmission
		map_indiv_stat = new HashMap<>();
		map_infhist_lines = new HashMap<>();

	}

//  For PID
//	int[] incl_start_grps = new int[] { 5, 6, 7, 8, 9 }; // Indigenous female
//	int[] sample_time = new int[] { 7300, 7665, 8030, 8395, 8760, 9125, 9490, 9855, 10220, 10585, 10950 };
//
//	int max_exposure = 120; // Assume won't develop PID after 4 months
//	double[] event_prob_by_inf_count = new double[] { 0.14, 0.17 };
//	int[] inf_count_range = new int[] { 0, 1 };

	public HashMap<String, double[]> analyse(int[] incl_start_grps, int[] sample_time, int max_exposure,
			double[] event_prob_by_inf_count, int[] inf_count_range) throws IOException {

		// Code start
		if (map_infhist_lines.isEmpty()) {
			loadInfHistMap(incl_start_grps);
		}
		HashMap<String, double[]> map_event_prob_by_file = new HashMap<>();

		for (Entry<String, ArrayList<int[]>> ent : map_infhist_lines.entrySet()) {
			Matcher m_map_infhist = infect_hist_key.matcher(ent.getKey());
			if (m_map_infhist.matches()) {
				Long cMap = Long.valueOf(m_map_infhist.group(3));
				HashMap<Integer, int[]> indivMap = map_indiv_stat.get(cMap);
				double[] event_prob_sum = new double[sample_time.length];
				map_event_prob_by_file.put(ent.getKey(), event_prob_sum);

				ArrayList<int[]> inf_hist_rows = ent.getValue();

				for (int r = 0; r < inf_hist_rows.size(); r++) {
					int[] row = inf_hist_rows.get(r);
					Integer id = row[0];
					int[] indiv_ent = indivMap.get(id);
					int start_grp = indiv_ent[Runnable_MetaPopulation_MultiTransmission.INDIV_MAP_ENTER_GRP];
					if (Arrays.binarySearch(incl_start_grps, start_grp) >= 0) {
						HashMap<Integer, ArrayList<Double>> map_cumul_prob_by_window_pt = new HashMap<>();

						int inf_count = 0;
						for (int inf_pt = 2; inf_pt < row.length; inf_pt += 3) { //
							int exposure_start = row[inf_pt];
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
								exposure_end = Math.min(exposure_end, row[inf_pt + 1]);
							}

							if (exposure_end > sample_time[0]) {
								int sample_end_point = Arrays.binarySearch(sample_time,
										Math.min(exposure_end, sample_time[sample_time.length - 1]));
								if (sample_end_point < 0) {
									sample_end_point = ~sample_end_point;
								}
								int sample_start_point = Arrays.binarySearch(sample_time,
										Math.max(sample_time[0], exposure_start));
								if (sample_start_point < 0) {
									sample_start_point = ~sample_start_point;
								}

								// # past infection
								int event_count_pt = Arrays.binarySearch(inf_count_range, inf_count);
								if (event_count_pt < 0) {
									event_count_pt = ~event_count_pt;
								}

								int current_window_pt = sample_start_point;

								while (current_window_pt <= sample_end_point) {
									if (current_window_pt > 0) {
										double current_event_prob = event_prob_by_inf_count[Math
												.min(event_prob_by_inf_count.length - 1, event_count_pt)];
										// P(0), P(1) etc... within current window
										ArrayList<Double> cumul_prob = map_cumul_prob_by_window_pt
												.get(current_window_pt);
										if (cumul_prob == null) {
											cumul_prob = new ArrayList<>();
											map_cumul_prob_by_window_pt.put(current_window_pt, cumul_prob);
											cumul_prob.add(Double.valueOf(1)); // P(0) = 1
										}

										int exposure_start_adj = Math.max(exposure_start,
												sample_time[current_window_pt - 1]);
										int exposure_end_adj = Math.min(exposure_end, sample_time[current_window_pt]);

										if (exposure_end_adj - exposure_start_adj < max_exposure) {
											double pDayExp = 1
													- Math.exp(Math.log(1 - current_event_prob) / max_exposure);
											current_event_prob = 1
													- Math.pow(1 - pDayExp, exposure_end_adj - exposure_start_adj);
										}

										// Update cumul_prod
										cumul_prob.add(cumul_prob.get(cumul_prob.size() - 1) * current_event_prob);
										for (int pI = cumul_prob.size() - 2; pI >= 1; pI--) {
											cumul_prob.set(pI, cumul_prob.get(pI) * (1 - current_event_prob)
													+ cumul_prob.get(pI - 1) * current_event_prob);
										}
										cumul_prob.set(0, cumul_prob.get(0) * (1 - current_event_prob));
									}
									current_window_pt++;

								}
							}
							inf_count++;
						}

						if (!map_cumul_prob_by_window_pt.isEmpty()) {
							Integer[] window_pts_array = map_cumul_prob_by_window_pt.keySet().toArray(new Integer[0]);
							Arrays.sort(window_pts_array);

							for (int wPt : window_pts_array) {
								ArrayList<Double> cumul_prob = map_cumul_prob_by_window_pt.get(wPt);
								if (cumul_prob.size() > 1) {
									// Update event_prob_sum
									double weighted_mean = 0;
									for (int i = 1; i < cumul_prob.size(); i++) {
										weighted_mean += i * cumul_prob.get(i);
									}
									event_prob_sum[wPt] += weighted_mean;
								}
							}
						}

					} // End of if (Arrays.binarySearch(incl_start_grps, start_grp) >= 0) {

				} // End of for (int r = 1; r < inf_hist_rows.size(); r++) {
			}
		}
		return map_event_prob_by_file;

	}

	private void loadInfHistMap(int[] incl_start_grps) throws IOException {

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
				// map_infhist_lines = StaticMethods.extractedLinesFrom7Zip(zip,
				// map_infhist_lines, null);
				SevenZArchiveEntry inputEnt;
				final int BUFFER = 2048;

				try {
					SevenZFile inputZip = new SevenZFile(zip);

					while ((inputEnt = inputZip.getNextEntry()) != null) {
						if (!inputEnt.isDirectory()) {
							String entName = inputEnt.getName();
							// Load map_indiv_stat for each map
							Matcher m_map_infhist = infect_hist_key.matcher(entName);
							if (m_map_infhist.matches()) {
								Long cMap = Long.valueOf(m_map_infhist.group(3));

								HashMap<Integer, int[]> indivMap = map_indiv_stat.get(cMap);
								if (indivMap == null) {
									indivMap = new HashMap<>();
									map_indiv_stat.put(cMap, indivMap);

									File demo_dir = new File(prop.getProperty("PROP_CONTACT_MAP_LOC"));
									demo_dir = new File(demo_dir, String.format("Demographic_%d", cMap));
									String[] pop_stat_line = StaticMethods.extracted_lines_from_text(
											new File(demo_dir, String.format("POP_STAT_%d.csv", cMap)));
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

								// Load all line
								int size = (int) inputEnt.getSize();

								byte[] content = new byte[size];
								int offset = 0;
								while (offset < size) {
									int readLen = inputZip.read(content, offset, Math.min(BUFFER, size - offset));
									if (readLen < 0) {
										break;
									}
									offset += readLen;
								}

								ArrayList<int[]> lines_split = new ArrayList<>();
								String[] lines = new String(content).split("\\n");

								for (int i = 1; i < lines.length; i++) {
									String line = lines[i];
									String[] lineEnt = line.split(",");
									int[] val = new int[lineEnt.length];
									for (int c = 0; c < val.length; c++) {
										val[c] = Integer.parseInt(lineEnt[c]);
									}
									if (incl_start_grps == null) {
										lines_split.add(val);
									} else {
										// Only load line with start_grp in incl_start_grps
										Integer id = val[0];
										int start_grp = indivMap
												.get(id)[Runnable_MetaPopulation_MultiTransmission.INDIV_MAP_ENTER_GRP];
										if (Arrays.binarySearch(incl_start_grps, start_grp) >= 0) {
											lines_split.add(val);
										}
									}

								}
								map_infhist_lines.put(entName, lines_split);

							} else {
								System.err.printf("Warning! Illformed zip file entry %s. Entry ignored.\n", entName);
							}

						}
					}
					inputZip.close();
				} catch (IOException e) {
					System.err.printf("Error when reading zip file %s. Zip file ignored.\n", zip.getAbsolutePath());
					e.printStackTrace(System.err);
				}

			}
		}

//		for (Entry<String, ArrayList<String[]>> ent : map_infhist_lines.entrySet()) {
//			Matcher m_map_infhist = infect_hist_key.matcher(ent.getKey());
//			if (m_map_infhist.matches()) {
//				Long cMap = Long.valueOf(m_map_infhist.group(3));
//
//				HashMap<Integer, int[]> indivMap = map_indiv_stat.get(cMap);
//				if (indivMap == null) {
//					indivMap = new HashMap<>();
//					map_indiv_stat.put(cMap, indivMap);
//
//					File demo_dir = new File(prop.getProperty("PROP_CONTACT_MAP_LOC"));
//					demo_dir = new File(demo_dir, String.format("Demographic_%d", cMap));
//					String[] pop_stat_line = StaticMethods
//							.extracted_lines_from_text(new File(demo_dir, String.format("POP_STAT_%d.csv", cMap)));
//					// ID,GRP,ENTER_POP_AGE,ENTER_POP_AT,EXIT_POP_AT,HOME_LOC
//					for (int i = 1; i < pop_stat_line.length; i++) {
//						String[] lineEnt = pop_stat_line[i].split(",");
//						int[] indiv_ent = new int[Runnable_MetaPopulation_MultiTransmission.LENGTH_INDIV_MAP];
//						indivMap.put(Integer.valueOf(lineEnt[0]), indiv_ent);
//						indiv_ent[Runnable_MetaPopulation_MultiTransmission.INDIV_MAP_ENTER_GRP] = Integer
//								.parseInt(lineEnt[1]);
//						indiv_ent[Runnable_MetaPopulation_MultiTransmission.INDIV_MAP_ENTER_POP_AGE] = Integer
//								.parseInt(lineEnt[2]);
//						indiv_ent[Runnable_MetaPopulation_MultiTransmission.INDIV_MAP_ENTER_POP_AT] = Integer
//								.parseInt(lineEnt[3]);
//						indiv_ent[Runnable_MetaPopulation_MultiTransmission.INDIV_MAP_EXIT_POP_AT] = Integer
//								.parseInt(lineEnt[4]);
//						indiv_ent[Runnable_MetaPopulation_MultiTransmission.INDIV_MAP_HOME_LOC] = Integer
//								.parseInt(lineEnt[5]);
//					}
//				}
//
//			} else {
//				System.err.printf("Warning! Illformed zip file entry %s. Entry ignored.\n", ent.getKey());
//			}
//		}
	}

	public static void generateInfectionHistCSV(HashMap<String, double[]> resultmap, int[] sample_time, File tarFile)
			throws FileNotFoundException {
		String[] resmap_key = resultmap.keySet().toArray(new String[0]);
		Arrays.sort(resmap_key, cmp_zipEnt);

		StringBuilder header = new StringBuilder();
		StringBuilder[] lines = new StringBuilder[sample_time.length];

		header.append("Time");
		for (int i = 0; i < lines.length; i++) {
			lines[i] = new StringBuilder();
			lines[i].append(sample_time[i]);
		}
		for (String zName : resmap_key) {
			Matcher m = patten_zipEnt.matcher(zName);
			m.matches();
			header.append(',');
			header.append(String.format("Seed_List_%s_%s", m.group(1), m.group(2)));
			double[] ent = resultmap.get(zName);
			for (int i = 0; i < ent.length; i++) {
				lines[i].append(',');
				lines[i].append(ent[i]);
			}

		}

		PrintWriter pWri = new PrintWriter(tarFile);
		pWri.println(header);
		for (StringBuilder line : lines) {
			pWri.println(line.toString());
		}

		pWri.close();
	}

}
