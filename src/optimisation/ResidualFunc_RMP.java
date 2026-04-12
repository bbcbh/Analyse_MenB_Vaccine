package optimisation;

import java.io.File;
import java.io.FileFilter;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.Properties;
import java.util.regex.Pattern;

import org.apache.commons.io.FileUtils;
import org.apache.commons.math3.analysis.MultivariateFunction;

import analysis.Analysis_PostSim_ExtractResults;
import sim.Abstract_Runnable_ClusterModel_Transmission;
import sim.Runnable_MenB_Vaccine;
import sim.Runnable_MenB_Vaccine_RMP;
import sim.Simulation_ClusterModelTransmission;
import sim.Simulation_MenB_Vaccine;

public class ResidualFunc_RMP implements MultivariateFunction {

	final String[] def_filepath;
	final String[] def_arg_simulation;
	final String[] def_arg_analysis;

	final String[] default_seed_file_header;
	final String[] default_seed_file_val;
	final String[] param_to_opt;
	final HashMap<String, Integer> opt_param_index_lookup;
	final HashMap<String, Integer> seed_list_param_index_lookup;

	final Map<String, String> param_cross_ref;

	final String[] opt_outcome_csv;
	// Format: double[OUT_COME_NUMBER]{target_val, conversion factor, weight}
	final double[][] opt_setting;
	final int[] opt_time_range;

	public static final int FILEPATH_SIM_DIR = 0;
	public static final int FILEPATH_SEED_DIR = 1;
	public static final int FILEPATH_REGION_MAP = 2;
	public static final int FILEPATH_GRP_SIZE = 3;

	public static final int DEFAULT_PARAMS_HEADER = 0;
	public static final int DEFAULT_PARAMS_INIT_VAL = 1;

	public static final int OPT_SETTING_TARGET = 0;
	public static final int OPT_SETTING_MODEL_POP_SIZE = 1;
	public static final int OPT_SETTING_WEIGHT = 2;

	public static final String fileformat_Opt_Outcomes = "OptProgress_ParamList_%s.csv";
	public static final String fileformat_Simplex_cache = "OptProgress_Simplex_%s.csv";
	public static final String fileformat_Output_txt = "Output.txt";

	protected final boolean printProgess;

	private HashMap<String, Double> eval_point_cache;
	private double minResidue = Double.POSITIVE_INFINITY;

	public ResidualFunc_RMP(String[] filepaths, String[][] default_params, String[] param_to_opt,
			Map<String, String> param_cross_ref, String[] opt_outcome_csv, double[][] opt_setting, int[] opt_time_range,
			boolean printProgess) {
		super();
		this.def_filepath = filepaths;
		this.default_seed_file_header = default_params[DEFAULT_PARAMS_HEADER];
		this.default_seed_file_val = default_params[DEFAULT_PARAMS_INIT_VAL];
		this.param_to_opt = param_to_opt;
		this.param_cross_ref = param_cross_ref;
		this.opt_outcome_csv = opt_outcome_csv;
		this.opt_setting = opt_setting;
		this.opt_time_range = opt_time_range;
		this.printProgess = printProgess;

		// Preset other fields
		this.def_arg_simulation = new String[] { filepaths[FILEPATH_SIM_DIR],
				"-seedMap=" + filepaths[FILEPATH_SEED_DIR] + File.separator + filepaths[FILEPATH_SEED_DIR] + ".csv",
				String.format("%s=%s", Simulation_ClusterModelTransmission.LAUNCH_ARGS_PRINT_PROGRESS, printProgess),
				"-export_skip_backup" };
		;
		this.def_arg_analysis = new String[] { filepaths[FILEPATH_SIM_DIR], filepaths[FILEPATH_REGION_MAP],
				filepaths[FILEPATH_GRP_SIZE],
				String.format("%s=%s", Simulation_ClusterModelTransmission.LAUNCH_ARGS_PRINT_PROGRESS, printProgess),
				"-flag=6" };

		// Use with seed_val_str
		seed_list_param_index_lookup = new HashMap<>();
		for (int i = 0; i < default_seed_file_header.length; i++) {
			seed_list_param_index_lookup.put(default_seed_file_header[i], i);
		}
		// Use with double[] point
		opt_param_index_lookup = new HashMap<>();
		for (int i = 0; i < param_to_opt.length; i++) {
			opt_param_index_lookup.put(param_to_opt[i], i);
		}

		// Write simplex values to cache
		eval_point_cache = new HashMap<>();
		File simplexCache = new File(def_filepath[FILEPATH_SIM_DIR],
				String.format(fileformat_Simplex_cache, def_filepath[FILEPATH_SEED_DIR]));
		if (simplexCache.exists()) {
			try {
				String[] ent = util.Util_7Z_CSV_Entry_Extract_Callable.extracted_lines_from_text(simplexCache);
				for (String line : ent) {
					String[] pair = line.split(":");
					eval_point_cache.put(pair[0], Double.valueOf(pair[1]));
				}
			} catch (IOException e) {
				e.printStackTrace(System.err);

			}
		}

	}

	@Override
	public double value(double[] point) {
		// Check if it is already cached
		if (eval_point_cache.containsKey(Arrays.toString(point))) {
			double res = eval_point_cache.get(Arrays.toString(point)).doubleValue();
			return res;
		} else {
			return value_eval(point);
		}
	}

	public double value_eval(double[] point) {
		double treatment_fit = Double.POSITIVE_INFINITY;

		String[] arg_simulation = Arrays.copyOf(def_arg_simulation, def_arg_simulation.length);
		String[] args_analysis = Arrays.copyOf(def_arg_analysis, def_arg_analysis.length);

		String[] seed_val_str = Arrays.copyOf(default_seed_file_val, default_seed_file_val.length);
		for (int i = 0; i < param_to_opt.length; i++) {
			String opt_param = param_to_opt[i];
			int seed_pt = seed_list_param_index_lookup.get(opt_param);
			double val = point[i];
			if (param_cross_ref.containsKey(opt_param)) {
				val *= point[opt_param_index_lookup.get(param_cross_ref.get(opt_param))];
			}
			seed_val_str[seed_pt] = Double.toString(val);
		}

		// Generate new working directory
		File working_dir = new File(def_filepath[FILEPATH_SIM_DIR], //
				String.format("OptDir_%s", def_filepath[FILEPATH_SEED_DIR]));
		if (working_dir.exists()) {
			// Remove previous working dir if already exist
			try {
				FileUtils.deleteDirectory(working_dir);
			} catch (Exception ex) {
				ex.printStackTrace(System.err);
			}
		}
		working_dir.mkdirs();

		File baseDir = new File(def_filepath[FILEPATH_SIM_DIR]);
		File[] copy_files = baseDir.listFiles(new FileFilter() {
			@Override
			public boolean accept(File pathname) {
				return pathname.isFile() && !pathname.getName().toLowerCase().startsWith("opt");
			}
		});

		arg_simulation[0] = working_dir.getAbsolutePath();
		args_analysis[0] = working_dir.getAbsolutePath();

		try {
			for (File f : copy_files) {
				Files.copy(f.toPath(), new File(working_dir, f.getName()).toPath(),
						StandardCopyOption.REPLACE_EXISTING);
			}
			// Generate seed file based on seed_val_str
			File seedDir = new File(working_dir, def_filepath[FILEPATH_SEED_DIR]);
			seedDir.mkdirs();
			PrintWriter pWri_seed = new PrintWriter(
					new File(seedDir, String.format("%s.csv", def_filepath[FILEPATH_SEED_DIR])));
			writeEntries(pWri_seed, default_seed_file_header);
			pWri_seed.println();
			writeEntries(pWri_seed, seed_val_str);
			pWri_seed.println();
			pWri_seed.close();

			// Run Simulation
			Simulation_MenB_Vaccine sim = new Simulation_MenB_Vaccine() {
				@Override
				public Abstract_Runnable_ClusterModel_Transmission generateDefaultRunnable(long cMap_seed,
						long sim_seed, Properties loadedProperties) {
					if (preGenSeedFile != null) {
						loadedProperties.put(Runnable_MenB_Vaccine.PROP_SEED_FILE_PATH,
								preGenSeedFile.getAbsolutePath());
					}
					Runnable_MenB_Vaccine_RMP runnable = new Runnable_MenB_Vaccine_RMP(cMap_seed, sim_seed,
							loadedProperties);
					return runnable;
				}

				@Override
				protected void finalise_simulations() throws IOException, FileNotFoundException {
					// Skip zipping of file
				}

			};
			sim.setPrintProgress(printProgess);
			Simulation_ClusterModelTransmission.launch(arg_simulation, sim);

			// Generate analyse result
			Analysis_PostSim_ExtractResults analyse_optRes = new Analysis_PostSim_ExtractResults();
			analyse_optRes.setMap_timetrend(Map.ofEntries( //
					Map.entry("Timetrend_Treatment_by_GrpLoc_%s.csv", Pattern.compile(".*Treatment_by_GrpLoc_.*csv"))));

			analyse_optRes.setRegion_extract_array(new String[][] { //
					new String[] { "Timetrend_Treatment_by_GrpLoc_%s.csv" }, //
					new String[] { "Timetrend_Treatment_by_GrpLoc_Grp_%d_Loc_%d.csv" }, //
					new String[] { "Timetrend_Treatment_%s_at_%s.csv" } });

			analyse_optRes.setPattern_sim(Pattern.compile(seedDir.getName()));

			analyse_optRes.analyse(args_analysis);

			// Calculate objective function based on analysis
			File outTxt = new File(working_dir, fileformat_Output_txt);
			PrintWriter wriOutTxt = new PrintWriter(new FileWriter(outTxt, true));

			wriOutTxt.printf("P=%s\n", Arrays.toString(seed_val_str));

			treatment_fit = 0;
			for (int f = 0; f < opt_outcome_csv.length; f++) {
				String[] lines = util.Util_7Z_CSV_Entry_Extract_Callable
						.extracted_lines_from_text(new File(working_dir, opt_outcome_csv[f]));
				int num_row_entries = 0;
				double residue_sum_by_outcome = 0;

				double treatment_rate_total = 0; // Might remove

				for (int r = 1; r < lines.length; r++) {
					double[] pre_rowEnt = null;
					String[] val = lines[r].split(",");
					int time = Integer.parseInt(val[0]);
					if (time >= opt_time_range[0] && time <= opt_time_range[1]) {
						if (pre_rowEnt == null) {
							String[] val_pre = lines[r - 1].split(",");
							pre_rowEnt = new double[val.length - 1];
							for (int c = 1; c < val.length; c++) {
								pre_rowEnt[c - 1] = Double.parseDouble(val_pre[c]);
							}
						}
						double[] rowEnt = new double[val.length - 1];
						for (int c = 1; c < val.length; c++) {
							rowEnt[c - 1] = Double.parseDouble(val[c]);
						}
						// Calculate sq diff
						double row_sq_diff = 0;
						for (int s = 0; s < rowEnt.length; s++) {
							double converted_rate = (rowEnt[s] - pre_rowEnt[s]) * 100000.0
									/ opt_setting[f][OPT_SETTING_MODEL_POP_SIZE];
							row_sq_diff += opt_setting[f][OPT_SETTING_WEIGHT]
									* Math.pow(opt_setting[f][OPT_SETTING_TARGET] - converted_rate, 2);
							treatment_rate_total += converted_rate;
						}

						row_sq_diff = row_sq_diff / rowEnt.length; // Average of all simulation (if more than one)
						residue_sum_by_outcome += row_sq_diff;

						pre_rowEnt = rowEnt;
						num_row_entries++;
					}
				}
				residue_sum_by_outcome = residue_sum_by_outcome / num_row_entries; // Average of included rows

				wriOutTxt.printf("#%d: Average treatment = %f from %d entries\n", f,
						treatment_rate_total / num_row_entries, num_row_entries);

				treatment_fit += residue_sum_by_outcome;
			}
			wriOutTxt.printf("R=%f\n", treatment_fit);			
			wriOutTxt.close();

			// Generate outcome file
			File file_outcome = new File(def_filepath[FILEPATH_SIM_DIR],
					String.format(fileformat_Opt_Outcomes, def_filepath[FILEPATH_SEED_DIR]));

			PrintWriter pWri_outcome;
			if (!file_outcome.exists()) {
				pWri_outcome = new PrintWriter(file_outcome);
				writeEntries(pWri_outcome, default_seed_file_header);
				pWri_outcome.println(",,RES");
			} else {
				pWri_outcome = new PrintWriter(new FileWriter(file_outcome, true));
			}
			writeEntries(pWri_outcome, seed_val_str);
			pWri_outcome.print(",,");
			pWri_outcome.print(treatment_fit);
			pWri_outcome.println();
			pWri_outcome.close();

			File file_simplexCache = new File(def_filepath[FILEPATH_SIM_DIR],
					String.format(fileformat_Simplex_cache, def_filepath[FILEPATH_SEED_DIR]));
			PrintWriter pWri_simplexCache = new PrintWriter(new FileWriter(file_simplexCache, true));
			pWri_simplexCache.printf("%s:%f\n", Arrays.toString(point), treatment_fit);
			pWri_simplexCache.close();

			if (treatment_fit < minResidue) {
				minResidue = treatment_fit;
				File target_dir = new File(working_dir.getParent(), String.format("BestFit_%s", working_dir.getName()));
				if (target_dir.exists()) {
					FileUtils.deleteDirectory(target_dir);
				}
				Files.move(working_dir.toPath(), target_dir.toPath());
			} else {
				// Remove working_dir recursively
				FileUtils.deleteDirectory(working_dir);
			}

		} catch (Exception e) {
			System.err.printf("Warning! %s encountered during running of parameter set: %s. Assume residue as Inf.\n",
					e.toString(), Arrays.toString(point));
			e.printStackTrace(System.err);
			System.exit(-1);
		}
		return treatment_fit;
	}

	private void writeEntries(PrintWriter pWri_seed, String[] arr) {
		for (int i = 0; i < arr.length; i++) {
			if (i != 0) {
				pWri_seed.append(',');
			}
			pWri_seed.append(arr[i]);
		}

	}

}
