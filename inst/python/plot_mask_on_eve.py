from nilearn import plotting, image
import nibabel as nib
import matplotlib.pyplot as plt
def plot_mask_on_eve(mask_img_path, template_img_path, cmap='bwr_r', alpha=1, save_path=None, **anat_kwargs):
    """
    Plot the mask image on the EVE template using Nilearn's plotting functionality.

    Parameters:
    - mask_img_path: str
        The file path to the mask NIfTI file.
    - template_img_path: str
        The file path to the template NIfTI file.
    - cmap: matplotlib.colors.Colormap, or str, optional
        The colormap to use. Either a string which is a name of a matplotlib colormap, 
        or a matplotlib colormap object.
    - alpha: float, optional
        The alpha transparency level for the overlay. Default is 1.
    - **anat_kwargs:
        Additional keyword arguments to pass to nilearn.plotting.plot_anat.
    """
    # Load the template and mask images
    template_img = nib.load(template_img_path)
    mask_img = nib.load(mask_img_path)

    # Adjust the affine matrix internally
    mask_img.affine[2][2] = 1
    mask_img.affine[0][0] = -1

    # Resample the mask to fit the template image
    mask_img_resampled = image.resample_to_img(mask_img, template_img, interpolation='nearest')

    # Plot the template and overlay the resampled mask
    display = plotting.plot_anat(template_img, **anat_kwargs)
    display.add_overlay(mask_img_resampled, cmap=cmap, alpha=alpha)

    # If a save path is provided, save the plot
    if save_path:
        plt.savefig(save_path)
        plt.close()  # Close the plot after saving to free up memory
    else:
        plt.show()  # Show the plot if not saving
    
    
