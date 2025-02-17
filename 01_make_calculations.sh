base=scripts

echo "======== Figure S10 ======="
python $base/figS10/calibration.py

echo "======== Figure 3 ======="
python $base/fig3/calc_visualization.py

echo "======== Figure S13 ======="
python $base/figS13/calc_visualization_threshold.py

echo "======== Figure S15 ======="
python $base/figS15/calc_noised_visualizations.py


