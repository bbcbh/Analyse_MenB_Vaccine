package optimisation;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.Properties;

import org.apache.commons.math3.analysis.MultivariateFunction;

public class MenB_RMP_NM_Optimistion extends Abstract_Optimisation {

	private final String path_region_path;
	private final String path_grp_size;

	private boolean objFunc_printProgress = false;

	public MenB_RMP_NM_Optimistion(String dirName, String seed_dir_name) throws IOException {
		super(dirName, seed_dir_name);

		// Additional Properties for MenB_RMP
		File file_opt_setting = new File(dirName, "optSetting.prop");
		FileInputStream fIS = new FileInputStream(file_opt_setting);
		Properties prop = new Properties();
		prop.loadFromXML(fIS);
		fIS.close();

		path_region_path = prop.getProperty("PROP_PATH_REGION_MAPPING");
		path_grp_size = prop.getProperty("PROP_PATH_GRP_SIZE");
		

	}

	protected MultivariateFunction generateObjectiveFunc(int seed_row, String[] seed_file_def_val) {
		String[] path = new String[] { simDirPath, String.format(OPTDIR_FORMAT, seedDirName, seed_row - 1),
				path_region_path, path_grp_size };
		ResidualFunc_RMP func = new ResidualFunc_RMP(path, new String[][] { seed_file_header, seed_file_def_val },
				param_to_opt, cross_ref_map, opt_outcome_csv, opt_setting, opt_time_range,
				objFunc_printProgress);
		return func;
	}

	public void setObjFunc_PrintProgress(boolean printProgress) {
		this.objFunc_printProgress = printProgress;
	}

}
