#!/usr/bin/env python3
"""
Priority Flow Workflow Example

This script walks through a typical workflow for processing a DEM to ensure a fully
connected hydrologic drainage network. For more details on the Priority Flow tool
refer to Condon and Maxwell (2019): https://doi.org/10.1016/j.cageo.2019.01.020

This is a Python translation of the R vignette Workflow_Example.Rmd from the
PriorityFlow R package (https://github.com/lecondon/PriorityFlow).

The example uses the sample watershed from Condon and Maxwell (2019). The DEM and
mask files are provided with the PriorityFlow library.
"""

# Plot inputs (optional - requires matplotlib)
def _plot_inputs(watershed_mask):
    """Plot the three input datasets."""
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    im0 = axes[0].imshow(DEM, cmap='RdBu' )
    axes[0].set_title("Elevation")
    plt.colorbar(im0, ax=axes[0])
    im1 = axes[1].imshow(watershed_mask, cmap='RdBu' )
    axes[1].set_title("Watershed Mask")
    plt.colorbar(im1, ax=axes[1])
    #im2 = axes[2].imshow(river_mask, cmap='RdBu' )
    #axes[2].set_title("River Network")
    #plt.colorbar(im2, ax=axes[2])
    plt.tight_layout()
    plt.savefig("workflow_inputs.png", dpi=150)
    plt.close()


def _plot_step1(trav_hs, dem_diff, targets, watershed_mask):
    """Plot DEM processing results."""
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    axes[0, 0].imshow(np.where(np.isnan(targets), 0, 1))
    axes[0, 0].set_title("Target Points")
    im1 = axes[0, 1].imshow(trav_hs["dem"])
    axes[0, 1].set_title("Processed DEM")
    plt.colorbar(im1, ax=axes[0, 1])
    im2 = axes[1, 0].imshow(dem_diff)
    axes[1, 0].set_title("DEM differences")
    plt.colorbar(im2, ax=axes[1, 0])
    im3 = axes[1, 1].imshow(trav_hs["direction"]/watershed_mask)
    axes[1, 1].set_title("Flow Direction")
    plt.colorbar(im3, ax=axes[1, 1])
    plt.tight_layout()
    plt.savefig("workflow_step1.png", dpi=150)
    plt.close()


def _plot_stream_network(
    subbasin: dict,
    stream_order: dict,
) -> None:
    """Plot stream segments, stream order, and subbasins (like R par(mfrow=c(1,3)))."""
    fig, axes = plt.subplots(1, 3, figsize=(14, 5))
    im0 = axes[0].imshow(subbasin["segments"], cmap="nipy_spectral")
    axes[0].set_title("Stream Segments")
    plt.colorbar(im0, ax=axes[0], label="Segment ID")
    im1 = axes[1].imshow(stream_order["order_mask"], cmap="viridis")
    axes[1].set_title("Stream Order")
    plt.colorbar(im1, ax=axes[1], label="Strahler order")
    im2 = axes[2].imshow(subbasin["subbasins"], cmap="nipy_spectral")
    axes[2].set_title("Subbasins")
    plt.colorbar(im2, ax=axes[2], label="Subbasin ID")
    plt.tight_layout()
    plt.savefig("workflow_stream_network.png", dpi=150)
    plt.close()

# Plot elevation differences from river smoothing
def _plot_step2():
    """Plot river smoothing results."""
    dif = riv_smooth_result["dem.adj"] - trav_hs["dem"]
    riv_mask = np.where(subbasin["segments"] > 0, 1, 0)
    hill_mask = 1 - riv_mask
    dif_hill = dif * hill_mask
    dif_riv = dif * riv_mask

    dif_plot = np.where(dif == 0, np.nan, dif)
    dif_riv_plot = np.where(dif_riv == 0, np.nan, dif_riv)
    dif_hill_plot = np.where(dif_hill == 0, np.nan, dif_hill)

    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    im0 = axes[0].imshow(dif_plot)
    axes[0].set_title("All Elev. Diffs")
    plt.colorbar(im0, ax=axes[0])
    im1 = axes[1].imshow(dif_riv_plot)
    axes[1].set_title("Stream Cell Diffs")
    plt.colorbar(im1, ax=axes[1])
    if np.any(~np.isnan(dif_hill_plot)):
        im2 = axes[2].imshow(dif_hill_plot)
        axes[2].set_title("Non-Stream Cell Diffs")
        plt.colorbar(im2, ax=axes[2])
    plt.tight_layout()
    plt.savefig("workflow_step2_smoothing.png", dpi=150)
    plt.close()


def _plot_path_transect(
    transect_old: np.ndarray,
    transect_new: np.ndarray,
    transect_riv: np.ndarray,
    subbasin: dict,
    streamline_riv: dict,
    segment: int,
) -> None:
    """
    Plot elevation transects and path map along a selected stream segment
    (Python analogue of the R PathExtract plotting block).
    """
    nstep = len(transect_riv)
    x = np.arange(1, nstep + 1)

    # Find breaks between stream segments along the path
    if nstep > 1:
        tr = np.asarray(transect_riv).ravel()
        diff_mask = tr[1:] != tr[:-1]
        slist = np.where(diff_mask)[0] + 1.5  # between steps k and k+1
    else:
        slist = []

    # Elevation limits across old/new transects
    all_vals = np.concatenate(
        [np.asarray(transect_old).ravel(), np.asarray(transect_new).ravel()]
    )
    vmin = np.nanmin(all_vals)
    vmax = np.nanmax(all_vals)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4))

    # Left: elevation transects
    ax0 = axes[0]
    ax0.plot(x, transect_old, color="C0", lw=2, label="Old Elevations")
    ax0.plot(x, transect_new, color="C2", lw=2, label="New Elevations")
    ax0.set_xlim(1, max(1, nstep))
    ax0.set_ylim(vmin, vmax)
    for s in slist:
        ax0.axvline(s, color="k", linestyle="--", linewidth=1)
    ax0.set_xlabel("Step")
    ax0.set_ylabel("Elevation")
    ax0.legend(loc="lower left", frameon=False)

    # Right: path map over stream network
    ax1 = axes[1]
    segment_plot = np.where(subbasin["segments"] > 0, 1, 0)
    ax1.imshow(segment_plot, cmap="gray_r", vmin=0, vmax=1)
    ax1.set_title(f"Path Map (segment {segment + 1})")

    path_mask = streamline_riv["path_mask"]
    path_overlay = np.where(path_mask == 0, np.nan, path_mask)
    ax1.imshow(path_overlay, cmap="Reds", alpha=0.8)

    plt.tight_layout()
    plt.savefig("workflow_path_transect.png", dpi=150)
    plt.close()


def _plot_slopes(slopex: np.ndarray, slopey: np.ndarray) -> None:
    """
    Plot resulting slopes in x and y directions (R analogue: image.plot of sxplot, syplot).
    """
    sxplot = np.where(slopex == 0, np.nan, slopex)
    syplot = np.where(slopey == 0, np.nan, slopey)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))

    im0 = axes[0].imshow(sxplot, cmap="RdBu")
    axes[0].set_title("SlopeX (primary & secondary)")
    plt.colorbar(im0, ax=axes[0])

    im1 = axes[1].imshow(syplot, cmap="RdBu")
    axes[1].set_title("SlopeY (primary & secondary)")
    plt.colorbar(im1, ax=axes[1])

    plt.tight_layout()
    plt.savefig("workflow_slopes.png", dpi=150)
    plt.close()

def write_parflow_ascii(data: np.ndarray, filepath: str) -> None:
    """
    Write a 2D array in ParFlow ASCII format.
    Format: header line with nx, ny, 1 followed by flattened data.
    """
    nx, ny = data.shape
    flat = np.zeros(nx * ny)
    jj = 0
    for j in range(ny):
        for i in range(nx):
            flat[jj] = data[i, j]
            jj += 1
    with open(filepath, "w") as f:
        f.write(f"{nx} {ny} 1\n")
        for val in flat:
            f.write(f"{val}\n")

