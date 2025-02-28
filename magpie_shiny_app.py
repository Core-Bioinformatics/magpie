import matplotlib.pyplot as plt
from shiny import App, ui, render, reactive, types
import numpy as np
from skimage.io import imread
import pandas as pd
import os
from shiny.types import ImgData
from pathlib import Path
from scipy import ndimage
import math
import glob
import copy

# Define the Shiny app UI
app_ui = ui.page_fluid(

    ui.tags.style(
        """
        .center-box {
            max-width: 800px;  /* Set the maximum width of the box */
            margin: 20px auto;    /* Center-aligns the box */
            background-color: #D1DAE2;
            padding: 15px;
            box-shadow: 0 4px 8px rgba(0,0,0,0.1); /* Optional shadow for styling */
            border-radius: 8px;
        }
        """
    ),
    
    ui.panel_well(
        
        ui.layout_columns(
            ui.panel_title("MAGPIE: Multimodal alignment of genes and peaks for integrative exploration"),
            ui.output_image("show_magpie_logo",height='300px'))),
    ui.panel_well(
        ui.h3("Sample selection:"),  
        ui.output_text("show_samples"),
        ui.input_select("pick_sample", None, choices=[],width='400px'),
        ui.output_text("show_peaks"),
        # Select dim reduction and run it
        ui.h3("Pick colouring for MSI data:"),
        ui.input_select("msi_colouring", None, choices=['PC1','First 3 PCs','Individual peak'],width='400px'),
        ui.panel_conditional("input.msi_colouring=='Individual peak'",
            ui.input_selectize("peak_choice", None, choices=[],multiple=False,width='400px')),
        ui.panel_conditional("input.msi_colouring!='Individual peak'",
            ui.input_selectize("peak_choices", None, choices=[],multiple=True,width='400px')),
        ui.input_slider("point_size","Select size of point",1,300,100,width='400px'),
        ui.h3("Adjust dimensionality reduction image"),
        ui.input_checkbox("flipx_dimred", "Flip dimensionality reduction in x", False),
        ui.input_checkbox("flipy_dimred", "Flip dimensionality reduction in y", False),
        ui.input_checkbox("rotate_dimred", "Rotate dimensionality reduction", False),
        ui.panel_conditional("input.rotate_dimred",
        ui.input_select("rotate_dimred_angle", "Angle to rotate dimensionality reduction", choices=[90,180,270])),
        ui.panel_conditional("output.middleHE=='MSI H&E image detected'",
            ui.h3("Adjust MSI H&E image"),
            ui.input_checkbox("flipx_msihe", "Flip MSI H&E in x", False),
            ui.input_checkbox("flipy_msihe", "Flip MSI H&E in y", False)),
            ui.input_checkbox("rotate_msihe", "Rotate MSI H&E", False),
            ui.panel_conditional("input.rotate_msihe",
            ui.input_select("rotate_msihe_angle", "Angle to rotate MSI H&E", choices=[90,180,270])),

        ui.input_action_button("run_dimred", "Run dimensionality reduction"),
        class_="center-box"),
    ui.panel_conditional("output.middleHE=='No MSI H&E image detected'",
                         ui.layout_columns(
                             ui.card(ui.output_plot("show_dim_red")),
                             ui.card(ui.output_plot("show_visium_he")),
                             row_heights="500px"
                         )),
    ui.panel_conditional("output.middleHE=='MSI H&E image detected'",
                        ui.layout_columns(
                            ui.card(ui.output_plot("show_dim_red2")),
                            ui.card(ui.output_plot("show_msi_he")),
                            ui.card(ui.output_plot("show_visium_he2")),
                            row_heights="500px"
                        )),                     
    # Show dim reduction
    ui.panel_conditional("input.flipx_dimred || input.flipy_dimred || input.rotate_dimred",
        ui.input_action_button("save_flipped_dim_red", "Save altered version of MSI data")),
    ui.panel_conditional("input.flipx_msihe || input.flipy_msihe || input.rotate_msihe",
        ui.input_action_button("save_flipped_msihe", "Save altered version of MSI H&E")),
    ui.input_action_button("save_dim_red", "Save MSI colouring for later analysis"),
    ui.output_text("middleHE"),
    # Pick landmarks between MSI and Visium H&E (if no MSI H&E)
    ui.panel_conditional("output.middleHE=='No MSI H&E image detected'",
                         ui.h3("Select landmarks between MSI data and Visium H&E"),
                         ui.layout_columns(
                            ui.card(ui.output_plot("plot_noHE_left", click=True)),  # Enable click on plot
                            ui.card(ui.output_plot("plot_noHE_right", click=True)),
                            row_heights="500px"),
                        ui.layout_columns(
                            ui.card(
                                ui.panel_conditional(
                                    "input.plot_noHE_left_click", 
                                    ui.output_plot("plot_noHE_withselected_left", click=False))),
                            ui.card(
                                ui.panel_conditional(
                                    "input.plot_noHE_right_click", 
                                    ui.output_plot("plot_noHE_withselected_right", click=False))),
                            row_heights="500px"),
                        ui.h4("Recorded Landmarks"),
                        ui.output_table("coords_noHE"),
                        ui.input_action_button("download_noHE", "Download landmarks"),
                        ui.input_action_button("undo_noHE_left_click", "Undo last MSI Image click"),
                        ui.input_action_button("undo_noHE_right_click", "Undo last Visium H&E click")),

    ui.panel_conditional("output.middleHE=='MSI H&E image detected'",
                        # Pick landmarks between MSI and MSI H&E
                        ui.h3("Select landmarks between MSI data and MSI H&E"),
                        ui.layout_columns(
                            ui.card(ui.output_plot("plot_MSI2HE_left", click=True)),  # Enable click on plot
                            ui.card(ui.output_plot("plot_MSI2HE_right", click=True)),
                            row_heights="500px"),
                        ui.layout_columns(
                            ui.card(
                                ui.panel_conditional(
                                    "input.plot_MSI2HE_left_click", 
                                    ui.output_plot("plot_MSI2HE_withselected_left", click=False))),
                            ui.card(
                                ui.panel_conditional(
                                    "input.plot_MSI2HE_right_click", 
                                    ui.output_plot("plot_MSI2HE_withselected_right", click=False))),
                            row_heights="500px"),
                        ui.h4("Recorded Landmarks"),
                        ui.output_table("coords_MSI2HE"),
                        ui.input_action_button("download_MSI2HE", "Download landmarks"),
                        ui.input_action_button("undo_MSI2HE_left_click", "Undo last MSI Image click"),
                        ui.input_action_button("undo_MSI2HE_right_click", "Undo last MSI H&E click"),
                        # Pick landmarks between MSI H&E and Visium H&E
                        ui.h2("Select landmarks between MSI H&E and Visium H&E"),
                         ui.layout_columns(
                            ui.card(ui.output_plot("plot_HE2HE_left", click=True)),  # Enable click on plot
                            ui.card(ui.output_plot("plot_HE2HE_right", click=True)),
                            row_heights="500px"),
                        ui.layout_columns(
                            ui.card(
                                ui.panel_conditional(
                                    "input.plot_HE2HE_left_click", 
                                    ui.output_plot("plot_HE2HE_withselected_left", click=False))),
                            ui.card(
                                ui.panel_conditional(
                                    "input.plot_HE2HE_right_click", 
                                    ui.output_plot("plot_HE2HE_withselected_right", click=False))),
                            row_heights="500px"),
                        ui.h4("Recorded Landmarks"),
                        ui.output_table("coords_HE2HE"),
                        ui.input_action_button("download_HE2HE", "Download landmarks"),
                        ui.input_action_button("undo_HE2HE_left_click", "Undo last MSI H&E click"),
                        ui.input_action_button("undo_HE2HE_right_click", "Undo last Visium H&E click")),
)

# Define the Shiny app server logic
def server(input, output, session):
    # Reactive value to store clicked coordinates
    clicked_coords_noHE_left = reactive.Value([])
    clicked_coords_noHE_right = reactive.Value([])
    clicked_coords_MSI2HE_left = reactive.Value([])
    clicked_coords_MSI2HE_right = reactive.Value([])
    clicked_coords_HE2HE_left = reactive.Value([])
    clicked_coords_HE2HE_right = reactive.Value([])

    @render.image
    def show_magpie_logo():
        dir = Path(__file__).resolve().parent
        img: ImgData = {"src": "figures/magpie_logo.png","height":"300px"}
        return img


    # Show available samples
    def find_samples():
        if os.path.isfile('input/selected.txt'):
            file = open("input/selected.txt", "r")
            visible_files = [line.rstrip() for line in file]
        else:
            visible_files = [x.replace('/msi', '').replace("\\msi","") for x in glob.glob('*/msi',root_dir='input')]
        if os.path.isfile('input/exclude.txt'):
            file = open("input/exclude.txt", "r")
            visible_files = list(set(visible_files).difference(set([line.rstrip() for line in file])))
        visible_files.sort()
        ui.update_select("pick_sample", choices=visible_files)
        return visible_files

    @output
    @render.text
    def show_samples():
        sample_list = find_samples()
        display_string = ", ".join(sample_list)
        display_string = "Found the following samples: " + display_string
        return(display_string)

    
    # Print whether or not there is an MSI H&E image (this seems to be necessary for the reactivity but could probably try to remove)
    @output
    @render.text
    @reactive.event(input.pick_sample)
    def middleHE():
        if glob.glob('input/'+input.pick_sample()+'/msi/MSI_HE.*') != []:
            msi_he_img = glob.glob('input/'+input.pick_sample()+'/msi/MSI_HE.*')
            msi_he_img = [x for x in msi_he_img if any(ext in x for ext in ['tiff','png','jpg'])]
            if not (msi_he_img is []):
                return('MSI H&E image detected')
            else:
                return('No MSI H&E image detected')
        else:
            return('No MSI H&E image detected')
    
    # Update choices for dim reduction based on whether there's an MSI H&E
    def dimred_options():
        msi_intensities = pd.read_csv('input/'+input.pick_sample()+'/msi/MSI_intensities.csv')
        msi_intensities.set_index('spot_id', drop=True, inplace=True)
        ui.update_selectize("peak_choice", choices=list(msi_intensities.columns))
        ui.update_selectize("peak_choices", choices=list(msi_intensities.columns))
        return list(msi_intensities.columns)
    
    @output
    @render.text
    def show_peaks():
        peak_list = dimred_options()
        display_string = "Found " + str(len(peak_list)) + " peaks in selected sample"
        return(display_string)

# Perform dimensionaly reduction
    @reactive.calc
    @reactive.event(input.run_dimred)
    def msi_dimred():

        clicked_coords_noHE_left.set([])
        clicked_coords_noHE_right.set([])
        clicked_coords_MSI2HE_left.set([])
        clicked_coords_MSI2HE_right.set([])
        clicked_coords_HE2HE_left.set([])
        clicked_coords_HE2HE_right.set([])
        
        msi_intensities = pd.read_csv('input/'+input.pick_sample()+'/msi/MSI_intensities.csv')
        msi_intensities.set_index('spot_id', drop=True, inplace=True)
        msi_coords = pd.read_csv('input/'+input.pick_sample()+'/msi/MSI_metadata.csv')
        if input.flipx_dimred()==True:
            msi_coords['x']= (-msi_coords['x'])
        if input.flipy_dimred()==True:
            msi_coords['y']= (-msi_coords['y'])
        if input.rotate_dimred()==True:

            old_x = msi_coords['x']
            old_y = msi_coords['y']
            msi_coords['x']= (old_x * math.cos(math.radians(int(input.rotate_dimred_angle())))) - (old_y * math.sin(math.radians(int(input.rotate_dimred_angle()))))
            msi_coords['y']= (old_x * math.sin(math.radians(int(input.rotate_dimred_angle())))) - (old_y * math.cos(math.radians(int(input.rotate_dimred_angle()))))

        if input.msi_colouring() == 'PC1':
            from sklearn.decomposition import PCA
            if list(input.peak_choices()) != []:
                msi_intensities = msi_intensities[list(input.peak_choices())]
            pca = PCA(n_components=1)
            reduction = pca.fit_transform(msi_intensities)
            msi_coords['color']=reduction
            
        if input.msi_colouring() == 'First 3 PCs':
            from sklearn.decomposition import PCA
            from sklearn.preprocessing import MinMaxScaler
            if list(input.peak_choices()) != []:
                msi_intensities = msi_intensities[list(input.peak_choices())]
            pca = PCA(n_components=3)
            reduction = pd.DataFrame(pca.fit_transform(msi_intensities))
            scaler = MinMaxScaler()
            reduction_scaled = pd.DataFrame(scaler.fit_transform(reduction), columns=reduction.columns)
            reduction_colours = reduction_scaled.values.tolist()
            def rgb_to_hex(r, g, b):
                return '#{:02x}{:02x}{:02x}'.format(r, g, b)
            reduction_colours_hex = [rgb_to_hex(int(np.round(255*x)),int(np.round(255*y)),int(np.round(255*z))) for [x,y,z] in reduction_colours]
            msi_coords['color']=reduction_colours_hex
#            ax.scatter(x=msi_coords['x'], y=msi_coords['y'], c=reduction_colours_hex,marker='.',s=input.point_size())
        if input.msi_colouring() == 'Individual peak':
            msi_coords['color'] = list(msi_intensities[input.peak_choice()])
        return(msi_coords)
    
    @reactive.Effect
    @reactive.event(input.save_dim_red)
    def save_flipped_dim_red():
        msi_coords = msi_dimred()
        msi_coords[['spot_id','x','y','color']].to_csv('input/'+input.pick_sample()+'/msi/MSI_dimreduction.csv',index=False)

    
    @reactive.calc
    @reactive.event(input.run_dimred)
    def msi_dimred_plot():
        msi_coords = msi_dimred()
        fig, ax = plt.subplots(nrows=1, ncols=1,dpi=100)  # create figure & 1 axis
        ax.margins(x=0,y=0)
        ax.scatter(x=msi_coords['x'], y=msi_coords['y'], c=msi_coords['color'],marker='.',s=input.point_size())
        fig.gca().set_aspect('equal')
        ax.set_title('MSI Image')
        ax.set_rasterization_zorder(0)
        fig.tight_layout()
        return (fig,ax)
    

    @output
    @render.plot(height=450)
    def show_dim_red():
        # Load the images
        fig,ax = copy.deepcopy(msi_dimred_plot())
        fig.set_dpi(100)
        return fig
    
    @output
    @render.plot(height=450)
    def show_dim_red2():
        # Load the images
        fig,ax = copy.deepcopy(msi_dimred_plot())
        fig.set_dpi(100)
        return fig
    
    # Show MSI H&E
    @output
    @render.plot(height=450)
    @reactive.event(input.pick_sample,input.flipx_msihe,input.flipy_msihe,input.rotate_msihe,input.rotate_msihe_angle)
    def show_msi_he():
        # Create the figure and axes
        fig, ax = plt.subplots(1, 1, figsize=(10, 5))
        # Load the images
        try:
            msi_he = imread(glob.glob('input/'+input.pick_sample()+'/msi/MSI_HE.*')[0])
            if input.flipx_msihe()==True:
                msi_he = np.fliplr(msi_he)
            if input.flipy_msihe()==True:
                msi_he = np.flipud(msi_he)
            if input.rotate_msihe()==True:
                msi_he = ndimage.rotate(msi_he, int(input.rotate_msihe_angle()))
            # Display the images in two subplots
            ax.imshow(msi_he)
            ax.set_title('MSI HE Image')

            # Tight layout for clean display
            fig.tight_layout()
            return fig
        except :
            pass
    
    @output
    @render.plot(height=450)
    @reactive.event(input.pick_sample)
    def show_visium_he():
        # Create the figure and axes
        fig, ax = plt.subplots(1, 1)
        # Load the images
        msi_he = imread('input/'+input.pick_sample()+'/visium/spatial/tissue_hires_image.png')
        
        # Display the images in two subplots
        ax.imshow(msi_he)
        ax.set_title('Visium HE Image')

        # Tight layout for clean display
        fig.tight_layout()
        return fig
    
    @output
    @render.plot(height=450)
    @reactive.event(input.pick_sample)
    def show_visium_he2():
        # Create the figure and axes
        fig, ax = plt.subplots(1, 1)
        # Load the images
        msi_he = imread('input/'+input.pick_sample()+'/visium/spatial/tissue_hires_image.png')
        
        # Display the images in two subplots
        ax.imshow(msi_he)
        ax.set_title('Visium HE Image')

        # Tight layout for clean display
        fig.tight_layout()
        return fig

    # All elements for picking landmarks between MSI and Visium H&E -------------------------

    # Show dim red and capture clicks
    @output
    @render.plot(height=450)
    def plot_noHE_left():
        # Load the images
        fig,ax = copy.deepcopy(msi_dimred_plot())
        fig.set_dpi(100)
        return fig
    
    # Show Visium H&E and capture clicks
    @output
    @render.plot(height=450)
    @reactive.event(input.pick_sample)
    def plot_noHE_right():
        # Create the figure and axes
        fig, ax = plt.subplots(1, 1)
        # Load the images
        msi_he = imread('input/'+input.pick_sample()+'/visium/spatial/tissue_hires_image.png')
        
        # Display the images in two subplots
        ax.imshow(msi_he)
        ax.set_title('Visium HE Image')

        # Tight layout for clean display
        fig.tight_layout()
        return fig

    # Show dim red with clicked landmarks
    @output
    @render.plot(height=450)
    @reactive.event(input.plot_noHE_left_click,input.undo_noHE_left_click)
    def plot_noHE_withselected_left():
        # Create the figure and axes
        fig,ax = copy.deepcopy(msi_dimred_plot())
#        # Get the list of clicked coordinates
        current_coords_left = clicked_coords_noHE_left.get()
        if current_coords_left:
            if len(current_coords_left)>0:
                x_vals, y_vals = zip(*current_coords_left)  # Unpack the coordinates
                for i in range(len(current_coords_left)):
                    ax.plot(x_vals[i], y_vals[i], 'ro', markersize=5)  # Red dots on MSI HE image
                    ax.text(x_vals[i], y_vals[i], str(i),   color='red',fontsize=9)
        fig.set_dpi(100)
        # Tight layout for clean display
        return fig
    
    # Show Visium H&E with clicked landmarks
    @output
    @render.plot(height=450)
    @reactive.event(input.plot_noHE_right_click,input.undo_noHE_right_click)
    def plot_noHE_withselected_right():
        # Create the figure and axes
        fig, ax = plt.subplots(1, 1, figsize=(10, 5))
        
        # Load the images
        msi_he = imread('input/'+input.pick_sample()+'/visium/spatial/tissue_hires_image.png')
       
        # Display the images in two subplots
        ax.imshow(msi_he)
        ax.set_title('Visium H&E Image')

        # Get the list of clicked coordinates
        current_coords_right = clicked_coords_noHE_right.get()
        if current_coords_right:
            if len(current_coords_right)>0:
                x_vals, y_vals = zip(*current_coords_right)  # Unpack the coordinates
                for i in range(len(current_coords_right)):
                    ax.plot(x_vals[i], y_vals[i], 'ro', markersize=5)  # Red dots on MSI HE image
                    ax.text(x_vals[i], y_vals[i], str(i),color='r',fontsize=9)
        # Tight layout for clean display
        fig.tight_layout()
        return fig

    # Update clicked points
    @reactive.Effect
    @reactive.event(input.plot_noHE_left_click)
    def update_noHE_leftclick():
        # Capture the click coordinates from the input
        click_info_left = input.plot_noHE_left_click()
        if click_info_left is not None:
            x = click_info_left['x']  # X-coordinate of the click
            y = click_info_left['y']  # Y-coordinate of the click

            # Update the list of clicked coordinates
            current_coords_left = clicked_coords_noHE_left.get()
            current_coords_left.append((x, y))  # Append the new coordinates
            clicked_coords_noHE_left.set(current_coords_left)

    # Undo clicked points
    @reactive.Effect
    @reactive.event(input.undo_noHE_left_click)
    def undo_noHE_leftclick():
        # Capture the click coordinates from the input
        current_coords_left = clicked_coords_noHE_left.get()
        current_coords_left = current_coords_left[:(len(current_coords_left)-1)]
        clicked_coords_noHE_left.set(current_coords_left)

    # Update clicked points
    @reactive.Effect
    @reactive.event(input.undo_noHE_right_click)
    def undo_noHE_rightclick():
        # Capture the click coordinates from the input
        current_coords_right = clicked_coords_noHE_right.get()
        current_coords_right = current_coords_right[:(len(current_coords_right)-1)]
        clicked_coords_noHE_right.set(current_coords_right)

    # Undo clicked points
    @reactive.Effect
    @reactive.event(input.plot_noHE_right_click)
    def update_noHE_rightclick():
        click_info_right = input.plot_noHE_right_click()
        if click_info_right is not None:
            x = click_info_right['x']  # X-coordinate of the click
            y = click_info_right['y']  # Y-coordinate of the click

            # Update the list of clicked coordinates
            current_coords_right = clicked_coords_noHE_right.get()
            current_coords_right.append((x, y))  # Append the new coordinates
            clicked_coords_noHE_right.set(current_coords_right)

    # Make coordinate table for landmarks MSI to Visium H&E
    @reactive.calc
    @reactive.event(input.plot_noHE_left_click,input.plot_noHE_right_click,input.undo_noHE_left_click,input.undo_noHE_right_click)
    def coords_calc_noHE():
        # Retrieve the stored coordinates
        current_coords_left = clicked_coords_noHE_left.get()
        current_coords_right = clicked_coords_noHE_right.get()
        
        if not current_coords_left:
            return pd.DataFrame(columns=["X_left", "Y_left", "X_right", "Y_right"])

        # Create a DataFrame to display the coordinates
        df_coords_left = pd.DataFrame(current_coords_left, columns=["X_left", "Y_left"])
        df_coords_right = pd.DataFrame(current_coords_right, columns=["X_right", "Y_right"])
        df_coords = pd.concat([df_coords_left,df_coords_right],axis=1)
        return df_coords
    
    # Display the clicked coordinates as a table
    @output
    @render.table
    def coords_noHE():
        return(coords_calc_noHE())

    @reactive.Effect
    @reactive.event(input.download_noHE)
    def download_noHE():
       coords_calc_noHE().to_csv('input/'+input.pick_sample()+'/landmarks_noHE.csv',index=False)

    @reactive.Effect
    @reactive.event(input.save_flipped_dim_red)
    def save_flipped_dim_red():
        msi_coords = pd.read_csv('input/'+input.pick_sample()+'/msi/MSI_metadata.csv')
        if input.flipx_dimred()==True:
            msi_coords['x']= (-msi_coords['x'])
        if input.flipy_dimred()==True:
            msi_coords['y']= (-msi_coords['y'])   
        if input.rotate_dimred()==True:
            old_x = msi_coords['x']
            old_y = msi_coords['y']
            msi_coords['x']= (old_x * math.cos(math.radians(input.rotate_dimred_angle()))) - (old_y * math.sin(math.radians(input.rotate_dimred_angle())))
            msi_coords['y']= (old_x * math.sin(math.radians(input.rotate_dimred_angle()))) - (old_y * math.cos(math.radians(input.rotate_dimred_angle())))

        msi_coords.to_csv('input/'+input.pick_sample()+'/msi/MSI_metadata_modified.csv',index=False)

    @reactive.Effect
    @reactive.event(input.save_flipped_msihe)
    def save_flipped_msihe():
        try:
            msi_he = imread(glob.glob('input/'+input.pick_sample()+'/msi/MSI_HE.*')[0])
            if input.flipx_msihe()==True:
                msi_he = np.fliplr(msi_he)
            if input.flipy_msihe()==True:
                msi_he = np.flipud(msi_he)
            if input.rotate_msihe()==True:
                msi_he = ndimage.rotate(msi_he, int(input.rotate_msihe_angle()))
            plt.imsave('input/'+input.pick_sample()+'/msi/MSI_HE_modified.jpg',msi_he)
        except:
            pass

     # All elements for picking landmarks between MSI and MSI H&E ---------------------------

    # Show dim red 
    @output
    @render.plot(height=450)
    def plot_MSI2HE_left():
        # Plot dim
        fig,ax = copy.deepcopy(msi_dimred_plot())
        fig.set_dpi(100)
        return fig
    
    # Show MSI H&E
    @output
    @render.plot(height=450)
    @reactive.event(input.pick_sample,input.flipx_msihe,input.flipy_msihe,input.rotate_msihe,input.rotate_msihe_angle)
    def plot_MSI2HE_right():
        # Create the figure and axes
        fig, ax = plt.subplots(1, 1, figsize=(10, 5))
        # Load the images
        try:
            msi_he = imread(glob.glob('input/'+input.pick_sample()+'/msi/MSI_HE.*')[0])
            if input.flipx_msihe()==True:
                msi_he = np.fliplr(msi_he)
            if input.flipy_msihe()==True:
                msi_he = np.flipud(msi_he)
            if input.rotate_msihe()==True:
                msi_he = ndimage.rotate(msi_he, int(input.rotate_msihe_angle()))
            # Display the images in two subplots
            ax.imshow(msi_he)
            ax.set_title('MSI HE Image')

            # Tight layout for clean display
            fig.tight_layout()
            return fig
        except:
            pass

    # Show dim red with clicked points
    @output
    @render.plot(height=450)
    @reactive.event(input.plot_MSI2HE_left_click,input.undo_MSI2HE_left_click)
    def plot_MSI2HE_withselected_left():
        # Create the figure and axes

        fig,ax = copy.deepcopy(msi_dimred_plot())

        # Get the list of clicked coordinates
        current_coords_left = clicked_coords_MSI2HE_left()
        if current_coords_left:
            if len(current_coords_left)>0:
                x_vals, y_vals = zip(*current_coords_left)  # Unpack the coordinates
                for i in range(len(current_coords_left)):
                    ax.plot(x_vals[i], y_vals[i], 'ro', markersize=5)  # Red dots on MSI HE image
                    ax.text(x_vals[i], y_vals[i], str(i),   color='red',fontsize=9)

        # Tight layout for clean display
        fig.set_dpi(100)
        return fig
    
    # Show MSI H&E with clicked points
    @output
    @render.plot(height=450)
    @reactive.event(input.plot_MSI2HE_right_click,input.undo_MSI2HE_right_click,input.flipx_msihe,input.flipy_msihe,input.rotate_msihe,input.rotate_msihe_angle)
    def plot_MSI2HE_withselected_right():
        # Create the figure and axes
        fig, ax = plt.subplots(1, 1, figsize=(10, 5))
        
        # Load the images
        try:
            msi_he = imread(glob.glob('input/'+input.pick_sample()+'/msi/MSI_HE.*')[0])
            if input.flipx_msihe()==True:
                msi_he = np.fliplr(msi_he)
            if input.flipy_msihe()==True:
                msi_he = np.flipud(msi_he)
            if input.rotate_msihe()==True:
                msi_he = ndimage.rotate(msi_he, int(input.rotate_msihe_angle()))

            # Display the images in two subplots
            ax.imshow(msi_he)
            ax.set_title('MSI H&E Image')

            # Get the list of clicked coordinates
            current_coords_right = clicked_coords_MSI2HE_right()
            if current_coords_right:
                if len(current_coords_right)>0:
                    x_vals, y_vals = zip(*current_coords_right)  # Unpack the coordinates
                    for i in range(len(current_coords_right)):
                        ax.plot(x_vals[i], y_vals[i], 'ro', markersize=5)  # Red dots on MSI HE image
                        ax.text(x_vals[i], y_vals[i], str(i),color='r',fontsize=9)
            # Tight layout for clean display
            fig.tight_layout()
            return fig
        except:
            pass

    # Update clicked points 
    @reactive.Effect
    @reactive.event(input.plot_MSI2HE_left_click)
    def update_MSI2HE_leftclick():
        # Capture the click coordinates from the input
        click_info_left = input.plot_MSI2HE_left_click()
        if click_info_left is not None:
            x = click_info_left['x']  # X-coordinate of the click
            y = click_info_left['y']  # Y-coordinate of the click

            # Update the list of clicked coordinates
            current_coords_left = clicked_coords_MSI2HE_left.get()
            current_coords_left.append((x, y))  # Append the new coordinates
            clicked_coords_MSI2HE_left.set(current_coords_left)

    # Update clicked points 
    @reactive.Effect
    @reactive.event(input.plot_MSI2HE_right_click)
    def update_MSI2HE_rightclick():
        click_info_right = input.plot_MSI2HE_right_click()
        if click_info_right is not None:
            x = click_info_right['x']  # X-coordinate of the click
            y = click_info_right['y']  # Y-coordinate of the click

            # Update the list of clicked coordinates
            current_coords_right = clicked_coords_MSI2HE_right.get()
            current_coords_right.append((x, y))  # Append the new coordinates
            clicked_coords_MSI2HE_right.set(current_coords_right)

    # Undo clicked points 
    @reactive.Effect
    @reactive.event(input.undo_MSI2HE_left_click)
    def undo_MSI2HE_leftclick():
        # Capture the click coordinates from the input
        current_coords_left = clicked_coords_MSI2HE_left.get()
        current_coords_left = current_coords_left[:(len(current_coords_left)-1)]
        clicked_coords_MSI2HE_left.set(current_coords_left)

    # Undo clicked points 
    @reactive.Effect
    @reactive.event(input.undo_MSI2HE_right_click)
    def undo_MSI2HE_rightclick():
        # Capture the click coordinates from the input
        current_coords_right = clicked_coords_MSI2HE_right.get()
        current_coords_right = current_coords_right[:(len(current_coords_right)-1)]
        clicked_coords_MSI2HE_right.set(current_coords_right)

    # Make table of landmarks
    @reactive.calc
    @reactive.event(input.plot_MSI2HE_left_click,input.plot_MSI2HE_right_click,input.undo_MSI2HE_left_click,input.undo_MSI2HE_right_click)
    def coords_calc_MSI2HE():
        # Retrieve the stored coordinates
        current_coords_left = clicked_coords_MSI2HE_left.get()
        current_coords_right = clicked_coords_MSI2HE_right.get()
        
        if not current_coords_left:
            return pd.DataFrame(columns=["X_left", "Y_left", "X_right", "Y_right"])

        # Create a DataFrame to display the coordinates
        df_coords_left = pd.DataFrame(current_coords_left, columns=["X_left", "Y_left"])
        df_coords_right = pd.DataFrame(current_coords_right, columns=["X_right", "Y_right"])
        df_coords = pd.concat([df_coords_left,df_coords_right],axis=1)
        return df_coords

    # Display the clicked coordinates as a table
    @output
    @render.table
    def coords_MSI2HE():
        # Retrieve the stored coordinates
        df_coords = coords_calc_MSI2HE()
        return df_coords
    
    # Download table
    @reactive.Effect
    @reactive.event(input.download_MSI2HE)
    def download_MSI2HE():
        coords_calc_MSI2HE().to_csv('input/'+input.pick_sample()+'/landmarks_MSI2HE.csv',index=False)


    # All elements for picking landmarks between MSI H&E and Visium ------------------------------------

    # Show MSI H&E
    @output
    @render.plot(height=450)
    @reactive.event(input.pick_sample,input.flipx_msihe,input.flipy_msihe,input.rotate_msihe,input.rotate_msihe_angle)
    def plot_HE2HE_left():
        # Load the images
       # Create the figure and axes
        fig, ax = plt.subplots(1, 1, figsize=(10, 5))
        # Load the images
        try:
            msi_he = imread(glob.glob('input/'+input.pick_sample()+'/msi/MSI_HE.*')[0])
            if input.flipx_msihe()==True:
                msi_he = np.fliplr(msi_he)
            if input.flipy_msihe()==True:
                msi_he = np.flipud(msi_he)
            if input.rotate_msihe()==True:
                msi_he = ndimage.rotate(msi_he, int(input.rotate_msihe_angle()))

            # Display the images in two subplots
            ax.imshow(msi_he)
            ax.set_title('MSI HE Image')

            # Tight layout for clean display
            fig.tight_layout()
            return fig
        except:
            pass
    
    # Show Visium H&E
    @output
    @render.plot(height=450)
    @reactive.event(input.pick_sample)
    def plot_HE2HE_right():
        # Create the figure and axes
        fig, ax = plt.subplots(1, 1, figsize=(10, 5))
        # Load the images
        msi_he = imread('input/'+input.pick_sample()+'/visium/spatial/tissue_hires_image.png')
        
        # Display the images in two subplots
        ax.imshow(msi_he)
        ax.set_title('Visium HE Image')

        # Tight layout for clean display
        fig.tight_layout()
        return fig

    # Show MSI H&E with clicked points
    @output
    @render.plot(height=450)
    @reactive.event(input.plot_HE2HE_left_click,input.undo_HE2HE_left_click,input.flipx_msihe,input.flipy_msihe,input.rotate_msihe,input.rotate_msihe_angle)
    def plot_HE2HE_withselected_left():
        # Create the figure and axes
        fig, ax = plt.subplots(1, 1, figsize=(10, 5))
        
        # Load the images
        try:
            msi_he = imread(glob.glob('input/'+input.pick_sample()+'/msi/MSI_HE.*')[0])
            if input.flipx_msihe()==True:
                msi_he = np.fliplr(msi_he)
            if input.flipy_msihe()==True:
                msi_he = np.flipud(msi_he)
            if input.rotate_msihe()==True:
                msi_he = ndimage.rotate(msi_he, int(input.rotate_msihe_angle()))

            # Display the images in two subplots
            ax.imshow(msi_he)
            ax.set_title('MSI H&E Image')

            # Get the list of clicked coordinates
            current_coords_left = clicked_coords_HE2HE_left()
            if current_coords_left:
                if len(current_coords_left)>0:
                    x_vals, y_vals = zip(*current_coords_left)  # Unpack the coordinates
                    for i in range(len(current_coords_left)):
                        ax.plot(x_vals[i], y_vals[i], 'ro', markersize=5)  # Red dots on MSI HE image
                        ax.text(x_vals[i], y_vals[i], str(i),color='r',fontsize=9)
            # Tight layout for clean display
            fig.tight_layout()
            return fig
        except:
            pass
    
    # Show Visium H&E with clicked points
    @output
    @render.plot(height=450)
    @reactive.event(input.plot_HE2HE_right_click,input.undo_HE2HE_right_click)
    def plot_HE2HE_withselected_right():
        # Create the figure and axes
        fig, ax = plt.subplots(1, 1, figsize=(10, 5))
        
        # Load the images
        visium_he = imread('input/'+input.pick_sample()+'/visium/spatial/tissue_hires_image.png')
       
        # Display the images in two subplots
        ax.imshow(visium_he)
        ax.set_title('Visium H&E Image')

        # Get the list of clicked coordinates
        current_coords_right = clicked_coords_HE2HE_right()
        if current_coords_right:
            if len(current_coords_right)>0:
                x_vals, y_vals = zip(*current_coords_right)  # Unpack the coordinates
                for i in range(len(current_coords_right)):
                    ax.plot(x_vals[i], y_vals[i], 'ro', markersize=5)  # Red dots on MSI HE image
                    ax.text(x_vals[i], y_vals[i], str(i),color='r',fontsize=9)
        # Tight layout for clean display
        fig.tight_layout()
        return fig

    # Update clicked points
    @reactive.Effect
    @reactive.event(input.plot_HE2HE_left_click)
    def update_HE2HE_leftclick():
        # Capture the click coordinates from the input
        click_info_left = input.plot_HE2HE_left_click()
        if click_info_left is not None:
            x = click_info_left['x']  # X-coordinate of the click
            y = click_info_left['y']  # Y-coordinate of the click

            # Update the list of clicked coordinates
            current_coords_left = clicked_coords_HE2HE_left.get()
            current_coords_left.append((x, y))  # Append the new coordinates
            clicked_coords_HE2HE_left.set(current_coords_left)

    # Update clicked points
    @reactive.Effect
    @reactive.event(input.plot_HE2HE_right_click)
    def update_HE2HE_rightclick():
        click_info_right = input.plot_HE2HE_right_click()
        if click_info_right is not None:
            x = click_info_right['x']  # X-coordinate of the click
            y = click_info_right['y']  # Y-coordinate of the click

            # Update the list of clicked coordinates
            current_coords_right = clicked_coords_HE2HE_right.get()
            current_coords_right.append((x, y))  # Append the new coordinates
            clicked_coords_HE2HE_right.set(current_coords_right)

    # Undo clicked points
    @reactive.Effect
    @reactive.event(input.undo_HE2HE_left_click)
    def undo_HE2HE_leftclick():
        # Capture the click coordinates from the input
        current_coords_left = clicked_coords_HE2HE_left.get()
        current_coords_left = current_coords_left[:(len(current_coords_left)-1)]
        clicked_coords_HE2HE_left.set(current_coords_left)

    # Undo clicked points
    @reactive.Effect
    @reactive.event(input.undo_HE2HE_right_click)
    def undo_HE2HE_rightclick():
        # Capture the click coordinates from the input
        current_coords_right = clicked_coords_HE2HE_right.get()
        current_coords_right = current_coords_right[:(len(current_coords_right)-1)]
        clicked_coords_HE2HE_right.set(current_coords_right)

    # Make table of landmarks
    @reactive.calc
    @reactive.event(input.plot_HE2HE_left_click,input.plot_HE2HE_right_click,input.undo_HE2HE_left_click,input.undo_HE2HE_right_click)
    def coords_calc_HE2HE():
        # Retrieve the stored coordinates
        current_coords_left = clicked_coords_HE2HE_left.get()
        current_coords_right = clicked_coords_HE2HE_right.get()
        
        if not current_coords_left:
            return pd.DataFrame(columns=["X_left", "Y_left", "X_right", "Y_right"])

        # Create a DataFrame to display the coordinates
        df_coords_left = pd.DataFrame(current_coords_left, columns=["X_left", "Y_left"])
        df_coords_right = pd.DataFrame(current_coords_right, columns=["X_right", "Y_right"])
        df_coords = pd.concat([df_coords_left,df_coords_right],axis=1)
        return df_coords
    
    # Display the clicked coordinates as a table
    @output
    @render.table
    def coords_HE2HE():
        # Retrieve the stored coordinates
        df_coords = coords_calc_HE2HE()
        return df_coords
        
    # Download table
    @reactive.Effect
    @reactive.event(input.download_HE2HE)
    def download_HE2HE():
        coords_calc_HE2HE().to_csv('input/'+input.pick_sample()+'/landmarks_HE2HE.csv',index=False)





# Create the Shiny app
app = App(app_ui, server)
