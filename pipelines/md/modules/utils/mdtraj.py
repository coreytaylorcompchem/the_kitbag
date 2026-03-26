import mdtraj as md

def reimage_trajectory(self, traj_path, top_path, output_path):

    logger.info(f"Re-imaging trajectory: {traj_path}")

    traj = md.load(traj_path, top=top_path)

    # Step 1: make molecules whole
    traj = traj.image_molecules()

    # Step 2: define centering selection
    top = traj.topology
    center_sel = top.select("protein or resname POP or resname UNK")

    if len(center_sel) == 0:
        logger.warning("Center selection empty - skipping centering")
    else:
        traj.center_coordinates(atom_indices=center_sel)

    # Step 3: re-image after centering
    traj = traj.image_molecules()

    traj.save(output_path)

    logger.info(f"Saved wrapped trajectory: {output_path}")