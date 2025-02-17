base=scripts

echo "======== Figure S10 ======="
python $base/figS10/figS10.py

echo "======== Figure 3 ======="
python $base/fig3/fig3a.py
python $base/fig3/fig3b.py
python $base/fig3/fig3c.py
python $base/fig3/fig3d.py
python $base/fig3/fig3ef.py

echo "======== Figure S12 ======="
python $base/figS12/figS12.py

echo "======== Figure S13 ======="
python $base/figS13/plot_visualizations_threshold.py

echo "======== Figure S14 ======="
python $base/figS14/plot_thresholds_interactions.py

echo "======== Figure S15 ======="
python $base/figS15/plot_noised_visualizations.py

echo "======== Figure S16 and S17 ======="
python $base/figS16_S17/figS16_S17.py

echo "======== Figure S18 ======="
python $base/figS18/figS18.py

