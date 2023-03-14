rm -fr ../output/*

alias ju="julia --project"

alias da="date \"+%Y-%m-%d %H:%M:%S\""

echo "⏳ $(da)"

# ---------------------------------------------------------------------------------------------- #

rename -e "s/_[\d]{4}-[\d]{2}-[\d]{2}_[\d]{2}-[\d]{2}-[\d]{2}//g" ../input/SDY*/* &&

ju 1.SDY1264.jl &&

ju 2.SDY67.jl &&

# ---------------------------------------------------------------------------------------------- #

echo "⌛️ $(da)"
