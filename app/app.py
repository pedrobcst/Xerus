import os
import shutil
from typing import List, Union

import pandas as pd
import streamlit as st
from st_aggrid import AgGrid
from Xerus.utils.cifutils import make_system_types

from conf import AppSettings
from xerus_plot import plot_highest_correlated, plot_read_data
from xerus_run import run_opt, run_xerus

os.makedirs(AppSettings.TMP_FOLDER, exist_ok=True)
os.makedirs(AppSettings.RESULTS_TMP_FOLDER, exist_ok=True)
from ui_utils import make_selection_label, process_input, read_input

st.title('XERUS Streamlit Interface Beta')
st.sidebar.markdown("**X**Ray **E**stimation and **R**efinement **U**sing **S**imilarity (**XERUS**)")
st.sidebar.image("https://raw.githubusercontent.com/pedrobcst/Xerus/master/img/g163.png", width=100)

base_columns = list(AppSettings.KEEP_COLUMNS)

# Session state stuff
if 'xerus_started' not in st.session_state:
    st.session_state['xerus_started'] = False

if 'xerus_object' not in st.session_state:
    st.session_state['xerus_object'] = None
if 'zip_file' not in st.session_state:
    st.session_state['zip_file'] = False

if 'optmized' not in st.session_state:
    st.session_state['optmized'] = False


@st.cache(allow_output_mutation=True)
def run_analysis(args_xerus: dict, args_analysis: dict):
    return run_xerus(args_xerus, args_analysis)

@st.cache(allow_output_mutation=True)
def run_optmizer(xerus_object, index_list: Union[int, List[int]], opt_args: dict):
    return run_opt(xerus_object, index_list, opt_args)


# Settings
st.sidebar.header("Data Upload and Name")
name = st.sidebar.text_input("Dataset name", key="name")
file = st.sidebar.file_uploader("Upload data", key="data_uploaded")
if file:
    data_format = st.sidebar.text_input("Data format", value=file.name.split(".")[-1], key="data_format")
    path = read_input(file)
    working_folder = os.path.join(AppSettings.RESULTS_TMP_FOLDER, file.name.split(".")[0]) + f"_{name}"
    os.makedirs(working_folder, exist_ok=True)

    # Data Visualization Settings
    with st.expander('Data View and Settings', expanded=True):
        if path:
            c1, c2 = st.columns([2, 4])
            # Background / Preprocessing information
            with c1:
                remove_background = st.checkbox("Remove background", value=True, key="remove_bg")
                use_preprocessed = st.checkbox("Use preprocessed data", value=False, key="use_pre")
                poly_degree = 10
                if remove_background:
                    poly_degree = st.number_input("Polynomial degree", min_value=2, max_value=12, step=1, value=8,key="poly_degree")    
                elements = st.text_input("Elements seperated by comma", value="Ho", key="element_list").split(",")
                elements = [element.strip() for element in elements if len(element) > 0]
                max_oxygen = st.number_input("Max oxygen", min_value=0, max_value=10, step=1, value=2, key="max_oxy", help="Define how much oxygen is allowed in the combinations. For example if max_oxy=2, and there are 3 elements (Ho,B,O), he will allow search up to Ho-O")
                is_ceramic = st.checkbox("Ceramic Material", key="is_ceramic", help = "If the material is a ceramic (oxide), the ignore_combs parameter later will be set to accept only oxide combinations.")
                st.write("Current element list is:", elements)
            # Data plot
            with c2:
                figure = plot_read_data(path, format=data_format, poly_degree=int(poly_degree), remove_base=remove_background)
                figure.update_layout(legend=dict(yanchor="top",y=1.35,xanchor="left",x=0.00, bgcolor="white"))
                st.plotly_chart(figure, use_container_width=True, template="presentation", bgcolor="white")
  
    c1, c2 = st.columns(2)
    # Algorithm settings
    with c1:
        with st.expander('Required Analysis Settings', expanded=True):
            n_runs = st.text_input("Number of runs", value="auto", key="n_runs")
            if n_runs != "auto":
                n_runs = int(n_runs)
            g = int(st.number_input("g", min_value=1, max_value=999, value=3, step=1, key="grabtop"))
            delta = st.number_input(r"delta", min_value=1.0, max_value=5.0, value=1.3, step=0.1, key="delta")
    # Optional search settings
    with c2:
        with st.expander("Optional Analysis Settings", expanded=True):
            ignore_ids = process_input(st.text_input("Ignore IDs", value="", key="ignore_ids"))
            ignore_providers = process_input(st.text_input("Ignore providers", value="AFLOW", key="ignore_providers"))
            ignore_comb = process_input(st.text_input("Ignore combinations", value="", key="ignore_comb"))
        if is_ceramic:
            to_ignore = make_system_types([ele for ele in elements if ele != "O"])
            if ignore_comb:
                ignore_comb.extend(to_ignore)
            else:
                ignore_comb = to_ignore
    
    # Filter Settings
    with st.expander("Current Filter Settings", expanded=False):
        c1, c2, c3 = st.columns(3)
        c1.markdown("**Ignore IDs**")
        c1.write(ignore_ids)
        c2.markdown("**Ignore combinations**")
        c2.write(ignore_comb)
        c3.markdown("**Ignore providers**")
        c3.write(ignore_providers)
    if st.button("Run analysis", key="run_analysis"):
        args_xerus = dict(name=name, working_folder=working_folder, exp_data_file=path, elements=elements,
                          max_oxy=max_oxygen, use_preprocessed=use_preprocessed, remove_background=remove_background,
                          poly_degree=poly_degree, data_fmt=data_format)

        args_analysis = dict(n_runs=n_runs, grabtop=g, delta=delta, ignore_ids=ignore_ids,
                             ignore_provider=ignore_providers, ignore_comb=ignore_comb)

        results_search = run_analysis(args_xerus, args_analysis)
        st.session_state['optmized'] = False
        st.session_state['zip_file'] = False
        st.balloons()
        st.session_state['xerus_object'] = results_search

    if st.session_state.xerus_object:
        st.header('Analysis Results')
        results_search = st.session_state.xerus_object
        df = results_search.results.copy()
        if 'wt' in df.columns:
            base_columns.append('wt')
        df = df[base_columns]
        df['id'] = df.index
        simuls_df = results_search.simulated_gsas2
        df = df[['id'] + base_columns]

        with st.expander('Results (Tabular)', expanded = True):
            AgGrid(df, width='50%', height=200)
        with st.expander('Visualization of Analysis', expanded = True):
            viz_number = st.selectbox("Results to Visualize", options=df.index, key="viz_number", format_func = lambda idx: make_selection_label(idx, df))
            fig = results_search.plot_result(viz_number)
            fig.update_layout(title=None, width=800, height=500)
            fig.update_xaxes(title=r'2theta (deg.)')
            st.plotly_chart(fig, use_container_width=True)
            plot_highest_corr = st.checkbox("Plot Highest correlated", value=False, key='plot_highest_corr')
            if plot_highest_corr:
                c1, c2 = st.columns([2,4])
                pattern_show = results_search.cif_info.copy()
                pattern_show.sort_values(by='Cij', inplace=True, ascending=False)
                pattern_show.reset_index(drop=True, inplace=True)

                highest_correlated = int(c1.number_input("Highest k correlated phases", min_value=1, max_value=len(simuls_df) - 1,
                                value=len(simuls_df) // 3, step=1, key='highest_corr'))

                options = pattern_show.filename[:highest_correlated]
                with c2:   
                    st.caption("List of all Options")
                    AgGrid(pattern_show.loc[:highest_correlated, ['filename', 'Cij', 'spacegroup']], width="50%", height=250)
                
                patterns_show = st.multiselect("Patterns to show", options=options, key='patterns_show', default=options)

                fig_highest_corr = plot_highest_correlated(data=results_search.exp_data_file, format=data_format,
                                                    cif_info=results_search.cif_info.copy(),
                                                    top=highest_correlated, width=800, height=500,
                                                    filter_columns=patterns_show)
                st.plotly_chart(fig_highest_corr, use_container_width=True)

        c1, c2 = st.columns(2)
        with c1:
            if st.button('Zip Contents'):
                shutil.make_archive(working_folder, 'zip', working_folder)
                st.session_state['zip_file'] = True
        with c2:
            if st.session_state['zip_file']:
                if os.path.exists(f"{working_folder}.zip"):
                    with open(f"{working_folder}.zip", "rb") as fp:
                        btn = st.download_button(
                        label="Download ZIP file",
                        data=fp,
                        file_name=f"{name}_results.zip",
                        mime="application/zip"
                    )
    
        with st.expander('Optimizer Settings'):
            # Basic optimizer settings
            c1,c2,c3,c4 = st.columns(4)
            optimizer_idx = process_input(c1.text_input('Indexes to optmize seperated by comma:', value="0", key="opt_list"),return_int=True)
            n_trials = int(c2.number_input("Number of trials", min_value=20, max_value=99999, value=200, step=1, key="n_trials"))
            param = c3.selectbox(label="Param to optimize", options=["rwp", "gof"])
            random_state = int(c4.number_input(label="Random seed number", min_value=0, max_value=9999, step=1, value=42,
                            key='random_state'))
            # Checkboxes
            c5, c6, c7, c8, c9 = st.columns(5)
            allow_pref_orient = c5.checkbox('Pref Orientation', value=True, key='pref_ori')
            allow_atomic_params = c6.checkbox('Atomic Params', value=False, key='atomic')
            allow_broad = c7.checkbox('Atomic Params', value=False, key='broadening')
            allow_angle = c8.checkbox('Acute angle', value=False, key='acute')
            force_ori = c9.checkbox('Force to use pref. ori', value=False, key='force_ori')

            opt_args = dict(n_trials=n_trials,
                            allow_pref_orient=allow_pref_orient,
                            allow_atomic_params=allow_atomic_params,
                            allow_broad=allow_broad,
                            allow_angle=allow_angle,
                            param=param,
                            random_state=random_state,
                            force_ori=force_ori
                            )
            st.write(opt_args)
            if st.button('Run optimization'):
                st.session_state['xerus_object'] = run_optmizer(results_search, optimizer_idx, opt_args)
                st.session_state['optmized'] = True
                st.balloons()
        if st.session_state['optmized']:
            with st.expander('Optimization Results'):
                st.write(
                 f'Optimization finished. Best rwp is {st.session_state.xerus_object.optimizer.optim.rwp_best:.3f}%')
                st.subheader('Refinement Result')
                fig = st.session_state.xerus_object.optimizer.optim.plot_best(save=False, engine="plotly")
                st.plotly_chart(fig, use_container_width=True)
                st.subheader('Crystal Structure Result')
                AgGrid(pd.DataFrame(data=st.session_state.xerus_object.optimizer.lattice_best), height=100)
                if st.button('Export Results to Working Folder'):
                    st.session_state.xerus_object.export_results()
                    st.write('Optimizaton results were exported to folder!')
                    st.write('Rezip and press download again!')

