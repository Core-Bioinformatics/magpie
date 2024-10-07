import matplotlib.pyplot as plt
from shiny import App, ui, render, reactive, types
import numpy as np
from skimage.io import imread
import pandas as pd
import os
from shiny.types import ImgData
from pathlib import Path


# Define the Shiny app UI
app_ui = ui.page_fluid(

    
    ui.panel_well(
        ui.layout_columns(
            ui.panel_title("MAGPIE: Multimodal alignment of genes and peaks for integrative exploration"),
            ui.output_image("show_magpie_logo",height='300px'))),
    
    ui.h3("Sample selection:"),  
    ui.output_text("show_samples"),
    ui.input_select("pick_sample", None, choices=[]),
    # Select dim reduction and run it
    ui.h3("Pick colouring for MSI data:"),
    ui.input_select("msi_colouring", None, choices=['PC1','First 3 PCs','Individual peak']),
    ui.panel_conditional("input.msi_colouring=='Individual peak'",
        ui.input_selectize("peak_choice", None, choices=[],multiple=False)),
    ui.panel_conditional("input.msi_colouring!='Individual peak'",
        ui.input_selectize("peak_choices", None, choices=[],multiple=True)),
    ui.input_slider("point_size","Select size of point",1,300,100),
    ui.input_action_button("run_dimred", "Run dimensionality reduction"),
    # Show dim reduction
    ui.output_plot("msi_dimred_output"),
    ui.output_text("middleHE"),
    # Pick landmarks between MSI and Visium H&E (if no MSI H&E)
    ui.panel_conditional("output.middleHE=='No MSI H&E image detected'",
                         ui.h3("Select landmarks between MSI data and Visium H&E"),
                         ui.layout_columns(
                            ui.card(ui.output_plot("plot_noHE_left", click=True),full_screen=True,height="400px"),  # Enable click on plot
                            ui.card(ui.output_plot("plot_noHE_right", click=True),full_screen=True,height="400px")),
                        ui.layout_columns(
                            ui.card(
                                ui.panel_conditional(
                                    "input.plot_noHE_left_click", 
                                    ui.output_plot("plot_noHE_withselected_left", click=False)),full_screen=True,height="400px"),
                            ui.card(
                                ui.panel_conditional(
                                    "input.plot_noHE_right_click", 
                                    ui.output_plot("plot_noHE_withselected_right", click=False)),full_screen=True,height="400px")),
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
                            ui.card(ui.output_plot("plot_MSI2HE_right", click=True))),
                        ui.layout_columns(
                            ui.card(
                                ui.panel_conditional(
                                    "input.plot_MSI2HE_left_click", 
                                    ui.output_plot("plot_MSI2HE_withselected_left", click=False))),
                            ui.card(
                                ui.panel_conditional(
                                    "input.plot_MSI2HE_right_click", 
                                    ui.output_plot("plot_MSI2HE_withselected_right", click=False)))),
                        ui.h4("Recorded Landmarks"),
                        ui.output_table("coords_MSI2HE"),
                        ui.input_action_button("download_MSI2HE", "Download landmarks"),
                        ui.input_action_button("undo_MSI2HE_left_click", "Undo last MSI Image click"),
                        ui.input_action_button("undo_MSI2HE_right_click", "Undo last MSI H&E click"),
                        # Pick landmarks between MSI H&E and Visium H&E
                        ui.h2("Select landmarks between MSI H&E and Visium H&E"),
                         ui.layout_columns(
                            ui.card(ui.output_plot("plot_HE2HE_left", click=True)),  # Enable click on plot
                            ui.card(ui.output_plot("plot_HE2HE_right", click=True))),
                        ui.layout_columns(
                            ui.card(
                                ui.panel_conditional(
                                    "input.plot_HE2HE_left_click", 
                                    ui.output_plot("plot_HE2HE_withselected_left", click=False))),
                            ui.card(
                                ui.panel_conditional(
                                    "input.plot_HE2HE_right_click", 
                                    ui.output_plot("plot_HE2HE_withselected_right", click=False)))),
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
        visible_files = sorted([file for file in os.listdir('input/') if not file.startswith('.')])
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
        if os.path.isfile('input/'+input.pick_sample()+'/msi/MSI_HE.jpg'):
            return('MSI H&E image detected')
        else :
            return('No MSI H&E image detected')
    
    # Update choices for dim reduction based on whether there's an MSI H&E
    @reactive.event(input.pick_sample)
    def dimred_options():
        msi_intensities = pd.read_csv('input/'+input.pick_sample()+'/msi/MSI_intensities.csv',index_col=0)
        ui.update_selectize("peak_choice", choices=list(msi_intensities.columns))
        ui.update_selectize("peak_choices", choices=list(msi_intensities.columns))
        return msi_intensities

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
        
        msi_intensities = dimred_options()
        msi_coords = pd.read_csv('input/'+input.pick_sample()+'/msi/MSI_metadata.csv',index_col=0)
        fig, ax = plt.subplots(nrows=1, ncols=1)  # create figure & 1 axis
        ax.margins(x=0,y=0)
        if input.msi_colouring() == 'PC1':
            from sklearn.decomposition import PCA
            if list(input.peak_choices()) != []:
                msi_intensities = msi_intensities[list(input.peak_choices())]
            pca = PCA(n_components=1)
            reduction = pca.fit_transform(msi_intensities)
            ax.scatter(x=msi_coords['x'], y=msi_coords['y'], c=reduction,marker='.',s=input.point_size())
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
            ax.scatter(x=msi_coords['x'], y=msi_coords['y'], c=reduction_colours_hex,marker='.',s=input.point_size())
        if input.msi_colouring() == 'Individual peak':
            ax.scatter(x=msi_coords['x'], y=msi_coords['y'], c=msi_intensities[input.peak_choice()],marker='.',s=input.point_size())
        fig.gca().set_aspect('equal')
        ax.set_title('MSI Image')
        ax.set_rasterization_zorder(0)
        plt.tight_layout()
        print('Im recalculating')
        return (fig,ax)

    # Show dimensionality reduction
    @output
    @render.plot
    def msi_dimred_output():
        fig,axs = msi_dimred()
        plt.tight_layout()
        return fig

    # All elements for picking landmarks between MSI and Visium H&E -------------------------

    # Show dim red and capture clicks
    @output
    @render.plot
    def plot_noHE_left():
        # Load the images
        fig,axs = msi_dimred()
        plt.tight_layout()
        return fig
    
    # Show Visium H&E and capture clicks
    @output
    @render.plot
    @reactive.event(input.pick_sample)
    def plot_noHE_right():
        # Create the figure and axes
        fig, axs = plt.subplots(1, 1)
        # Load the images
        msi_he = imread('input/'+input.pick_sample()+'/visium/spatial/tissue_hires_image.png')
        
        # Display the images in two subplots
        axs.imshow(msi_he)
        axs.set_title('Visium HE Image')

        # Tight layout for clean display
        plt.tight_layout()
        return fig

    # Show dim red with clicked landmarks
    @output
    @render.plot
    @reactive.event(input.plot_noHE_left_click,input.undo_noHE_left_click)
    def plot_noHE_withselected_left():
        # Create the figure and axes
        fig, axs = msi_dimred()
#        # Get the list of clicked coordinates
        current_coords_left = clicked_coords_noHE_left.get()
        if current_coords_left:
            if len(current_coords_left)>0:
                x_vals, y_vals = zip(*current_coords_left)  # Unpack the coordinates
                for i in range(len(current_coords_left)):
                    axs.plot(x_vals[i], y_vals[i], 'ro', markersize=5)  # Red dots on MSI HE image
                    axs.text(x_vals[i], y_vals[i], str(i),   color='red',fontsize=9)
        # Tight layout for clean display
        return fig
    
    # Show Visium H&E with clicked landmarks
    @output
    @render.plot
    @reactive.event(input.plot_noHE_right_click,input.undo_noHE_right_click)
    def plot_noHE_withselected_right():
        # Create the figure and axes
        fig, axs = plt.subplots(1, 1, figsize=(10, 5))
        
        # Load the images
        msi_he = imread('input/'+input.pick_sample()+'/visium/spatial/tissue_hires_image.png')
       
        # Display the images in two subplots
        axs.imshow(msi_he)
        axs.set_title('Visium H&E Image')

        # Get the list of clicked coordinates
        current_coords_right = clicked_coords_noHE_right.get()
        if current_coords_right:
            if len(current_coords_right)>0:
                x_vals, y_vals = zip(*current_coords_right)  # Unpack the coordinates
                for i in range(len(current_coords_right)):
                    axs.plot(x_vals[i], y_vals[i], 'ro', markersize=5)  # Red dots on MSI HE image
                    axs.text(x_vals[i], y_vals[i], str(i),color='r',fontsize=9)
        # Tight layout for clean display
        plt.tight_layout()
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
       coords_calc_noHE().to_csv('input/'+input.pick_sample()+'/landmarks_noHE.csv')

     # All elements for picking landmarks between MSI and MSI H&E ---------------------------

    # Show dim red 
    @output
    @render.plot
    def plot_MSI2HE_left():
        # Plot dim
        fig,axs = msi_dimred()
        return fig
    
    # Show MSI H&E
    @output
    @render.plot
    @reactive.event(input.pick_sample)
    def plot_MSI2HE_right():
        # Create the figure and axes
        fig, axs = plt.subplots(1, 1, figsize=(10, 5))
        # Load the images
        msi_he = imread('input/'+input.pick_sample()+'/msi/MSI_HE.jpg')
        
        # Display the images in two subplots
        axs.imshow(msi_he)
        axs.set_title('MSI HE Image')

        # Tight layout for clean display
        plt.tight_layout()
        return fig

    # Show dim red with clicked points
    @output
    @render.plot
    @reactive.event(input.plot_MSI2HE_left_click,input.undo_MSI2HE_left_click)
    def plot_MSI2HE_withselected_left():
        # Create the figure and axes

        fig, axs = msi_dimred()

        # Get the list of clicked coordinates
        current_coords_left = clicked_coords_MSI2HE_left()
        if current_coords_left:
            if len(current_coords_left)>0:
                x_vals, y_vals = zip(*current_coords_left)  # Unpack the coordinates
                for i in range(len(current_coords_left)):
                    axs.plot(x_vals[i], y_vals[i], 'ro', markersize=5)  # Red dots on MSI HE image
                    axs.text(x_vals[i], y_vals[i], str(i),   color='red',fontsize=9)

        # Tight layout for clean display
        fig.tight_layout()
        return fig
    
    # Show MSI H&E with clicked points
    @output
    @render.plot
    @reactive.event(input.plot_MSI2HE_right_click,input.undo_MSI2HE_right_click)
    def plot_MSI2HE_withselected_right():
        # Create the figure and axes
        fig, axs = plt.subplots(1, 1, figsize=(10, 5))
        
        # Load the images

        msi_he = imread('input/'+input.pick_sample()+'/msi/MSI_HE.jpg')
       
        # Display the images in two subplots
        axs.imshow(msi_he)
        axs.set_title('MSI H&E Image')

        # Get the list of clicked coordinates
        current_coords_right = clicked_coords_MSI2HE_right()
        if current_coords_right:
            if len(current_coords_right)>0:
                x_vals, y_vals = zip(*current_coords_right)  # Unpack the coordinates
                for i in range(len(current_coords_right)):
                    axs.plot(x_vals[i], y_vals[i], 'ro', markersize=5)  # Red dots on MSI HE image
                    axs.text(x_vals[i], y_vals[i], str(i),color='r',fontsize=9)
        # Tight layout for clean display
        plt.tight_layout()
        return fig

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
        coords_calc_MSI2HE().to_csv('input/'+input.pick_sample()+'/landmarks_MSI2HE.csv')


    # All elements for picking landmarks between MSI H&E and Visium ------------------------------------

    # Show MSI H&E
    @output
    @render.plot
    @reactive.event(input.pick_sample)
    def plot_HE2HE_left():
        # Load the images
       # Create the figure and axes
        fig, axs = plt.subplots(1, 1, figsize=(10, 5))
        # Load the images
        msi_he = imread('input/'+input.pick_sample()+'/msi/MSI_HE.jpg')
        
        # Display the images in two subplots
        axs.imshow(msi_he)
        axs.set_title('MSI HE Image')

        # Tight layout for clean display
        plt.tight_layout()
        return fig
    
    # Show Visium H&E
    @output
    @render.plot
    @reactive.event(input.pick_sample)
    def plot_HE2HE_right():
        # Create the figure and axes
        fig, axs = plt.subplots(1, 1, figsize=(10, 5))
        # Load the images
        msi_he = imread('input/'+input.pick_sample()+'/visium/spatial/tissue_hires_image.png')
        
        # Display the images in two subplots
        axs.imshow(msi_he)
        axs.set_title('Visium HE Image')

        # Tight layout for clean display
        plt.tight_layout()
        return fig

    # Show MSI H&E with clicked points
    @output
    @render.plot
    @reactive.event(input.plot_HE2HE_left_click,input.undo_HE2HE_left_click)
    def plot_HE2HE_withselected_left():
        # Create the figure and axes
        fig, axs = plt.subplots(1, 1, figsize=(10, 5))
        
        # Load the images
        msi_he = imread('input/'+input.pick_sample()+'/msi/MSI_HE.jpg')
       
        # Display the images in two subplots
        axs.imshow(msi_he)
        axs.set_title('MSI H&E Image')

        # Get the list of clicked coordinates
        current_coords_left = clicked_coords_HE2HE_left()
        if current_coords_left:
            if len(current_coords_left)>0:
                x_vals, y_vals = zip(*current_coords_left)  # Unpack the coordinates
                for i in range(len(current_coords_left)):
                    axs.plot(x_vals[i], y_vals[i], 'ro', markersize=5)  # Red dots on MSI HE image
                    axs.text(x_vals[i], y_vals[i], str(i),color='r',fontsize=9)
        # Tight layout for clean display
        plt.tight_layout()
        return fig
    
    # Show Visium H&E with clicked points
    @output
    @render.plot
    @reactive.event(input.plot_HE2HE_right_click,input.undo_HE2HE_right_click)
    def plot_HE2HE_withselected_right():
        # Create the figure and axes
        fig, axs = plt.subplots(1, 1, figsize=(10, 5))
        
        # Load the images
        visium_he = imread('input/'+input.pick_sample()+'/visium/spatial/tissue_hires_image.png')
       
        # Display the images in two subplots
        axs.imshow(visium_he)
        axs.set_title('Visium H&E Image')

        # Get the list of clicked coordinates
        current_coords_right = clicked_coords_HE2HE_right()
        if current_coords_right:
            if len(current_coords_right)>0:
                x_vals, y_vals = zip(*current_coords_right)  # Unpack the coordinates
                for i in range(len(current_coords_right)):
                    axs.plot(x_vals[i], y_vals[i], 'ro', markersize=5)  # Red dots on MSI HE image
                    axs.text(x_vals[i], y_vals[i], str(i),color='r',fontsize=9)
        # Tight layout for clean display
        plt.tight_layout()
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
        coords_calc_HE2HE().to_csv('input/'+input.pick_sample()+'/landmarks_HE2HE.csv')





# Create the Shiny app
app = App(app_ui, server)
