%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%	MDPI class for LaTeX files	 03-03-2015
%% %%	For any information please send an e-mail to:
%% %%		latex-support@mdpi.com (Janine Daum) 
%% %%
%% %%	Class style provided by Stefano Mariani (stefano.mariani@polimi.it)
%% %%   Modified by:
%% %%		Dietrich Rordorf (rordorf@mdpi.com)
%% %%		Peter Harremoes (entropeter@mdpi.com)
%% %%		Zeno Schumacher (latex@mdpi.com)
%% %%		Maddalena Giulini (latex@mdpi.com)
%% %%		Andres Gartmann (gartmann@mdpi.com) v1.2
%% %%		Janine Daum (latex-support@mdpi.com) v3.0
%% %%		
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% IDENTIFICATION
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{mdpi}[2015/03/03 MDPI paper class]

\AtEndDocument{\par \cright} %% needs to be here

%% PRELIMINARY DECLARATIONS
\LoadClass[10pt]{article}
\RequirePackage[T1]{fontenc}
\RequirePackage[utf8]{inputenc}
\RequirePackage{calc}
\RequirePackage{indentfirst}
\RequirePackage{fancyhdr}
\RequirePackage{graphicx,epstopdf}
\RequirePackage{lastpage}
\RequirePackage{ifthen}
\RequirePackage{lineno}
\RequirePackage{float}
\RequirePackage{amsmath}
\RequirePackage{setspace}
\RequirePackage{enumitem}
\RequirePackage{mathpazo}
\RequirePackage{booktabs} % for \toprule etc. in tables
\RequirePackage[largestsep]{titlesec}
\RequirePackage{etoolbox} % for \AtBeginDocument

%% OPTIONS
%% To choose the journal
\DeclareOption{actuators}{\gdef\@journal{actuators}}
\DeclareOption{admsci}{\gdef\@journal{admsci}}
\DeclareOption{administrativesciences}{\gdef\@journal{admsci}} %to have backwards compatibility
\DeclareOption{aerospace}{\gdef\@journal{aerospace}}
\DeclareOption{agriculture}{\gdef\@journal{agriculture}}
\DeclareOption{agronomy}{\gdef\@journal{agronomy}}
\DeclareOption{algorithms}{\gdef\@journal{algorithms}}
\DeclareOption{animals}{\gdef\@journal{animals}}
\DeclareOption{antibiotics}{\gdef\@journal{antibiotics}}
\DeclareOption{antibodies}{\gdef\@journal{antibodies}}
\DeclareOption{antioxidants}{\gdef\@journal{antioxidants}}
\DeclareOption{applsci}{\gdef\@journal{applsci}}
\DeclareOption{appliedsciences}{\gdef\@journal{applsci}} %to have backwards compatibility
\DeclareOption{arts}{\gdef\@journal{arts}}
\DeclareOption{atmosphere}{\gdef\@journal{atmosphere}}
\DeclareOption{atoms}{\gdef\@journal{atoms}}
\DeclareOption{axioms}{\gdef\@journal{axioms}}
\DeclareOption{batteries}{\gdef\@journal{batteries}}
\DeclareOption{behavsci}{\gdef\@journal{behavsci}} 
\DeclareOption{behavioralsciences}{\gdef\@journal{behavsci}} %to have backwards compatibility 
\DeclareOption{beverages}{\gdef\@journal{beverages}}
\DeclareOption{bioengineering}{\gdef\@journal{bioengineering}}
\DeclareOption{biology}{\gdef\@journal{biology}}
\DeclareOption{biomedicines}{\gdef\@journal{biomedicines}}
\DeclareOption{biomimetics}{\gdef\@journal{biomimetics}}
\DeclareOption{biomolecules}{\gdef\@journal{biomolecules}}
\DeclareOption{biosensors}{\gdef\@journal{biosensors}}
\DeclareOption{brainsci}{\gdef\@journal{brainsci}} 
\DeclareOption{brainsciences}{\gdef\@journal{brainsci}} %to have backwards compatibility 
\DeclareOption{buildings}{\gdef\@journal{buildings}}
\DeclareOption{c}{\gdef\@journal{carbon}}%% for authors
\DeclareOption{carbon}{\gdef\@journal{carbon}} 
\DeclareOption{cancers}{\gdef\@journal{cancers}}
\DeclareOption{catalysts}{\gdef\@journal{catalysts}}
\DeclareOption{cells}{\gdef\@journal{cells}}
\DeclareOption{challenges}{\gdef\@journal{challenges}}
\DeclareOption{chemosensors}{\gdef\@journal{chemosensors}}
\DeclareOption{children}{\gdef\@journal{children}}
\DeclareOption{chromatography}{\gdef\@journal{chromatography}}
\DeclareOption{climate}{\gdef\@journal{climate}}
\DeclareOption{coatings}{\gdef\@journal{coatings}}
\DeclareOption{computation}{\gdef\@journal{computation}}
\DeclareOption{computers}{\gdef\@journal{computers}}
\DeclareOption{condensedmatter}{\gdef\@journal{condensedmatter}}
\DeclareOption{cosmetics}{\gdef\@journal{cosmetics}}
\DeclareOption{cryptography}{\gdef\@journal{cryptography}}
\DeclareOption{crystals}{\gdef\@journal{crystals}}
\DeclareOption{data}{\gdef\@journal{data}}
\DeclareOption{dentistry}{\gdef\@journal{dentistry}}
\DeclareOption{dentistryjournal}{\gdef\@journal{dentistry}} %to have backwards compatibility 
\DeclareOption{designs}{\gdef\@journal{designs}} 
\DeclareOption{diagnostics}{\gdef\@journal{diagnostics}}
\DeclareOption{diseases}{\gdef\@journal{diseases}}
\DeclareOption{diversity}{\gdef\@journal{diversity}}
\DeclareOption{econometrics}{\gdef\@journal{econometrics}}
\DeclareOption{economies}{\gdef\@journal{economies}}
\DeclareOption{education}{\gdef\@journal{education}}
\DeclareOption{educationsciences}{\gdef\@journal{education}} %to have backwards compatibility 
\DeclareOption{electronics}{\gdef\@journal{electronics}}
\DeclareOption{energies}{\gdef\@journal{energies}}
\DeclareOption{entropy}{\gdef\@journal{entropy}}
\DeclareOption{environments}{\gdef\@journal{environments}}
\DeclareOption{epigenomes}{\gdef\@journal{epigenomes}}
\DeclareOption{fermentation}{\gdef\@journal{fermentation}}
\DeclareOption{fibers}{\gdef\@journal{fibers}}
\DeclareOption{fishes}{\gdef\@journal{fishes}}
\DeclareOption{fluids}{\gdef\@journal{fluids}}
\DeclareOption{foods}{\gdef\@journal{foods}}
\DeclareOption{forests}{\gdef\@journal{forests}}
\DeclareOption{futureinternet}{\gdef\@journal{futureinternet}}
\DeclareOption{galaxies}{\gdef\@journal{galaxies}}
\DeclareOption{games}{\gdef\@journal{games}}
\DeclareOption{gels}{\gdef\@journal{gels}}
\DeclareOption{genealogy}{\gdef\@journal{genealogy}}
\DeclareOption{genes}{\gdef\@journal{genes}}
\DeclareOption{geosciences}{\gdef\@journal{geosciences}}
\DeclareOption{geriatrics}{\gdef\@journal{geriatrics}}
\DeclareOption{healthcare}{\gdef\@journal{healthcare}}
\DeclareOption{horticulturae}{\gdef\@journal{horticulturae}}
\DeclareOption{humanities}{\gdef\@journal{humanities}}
\DeclareOption{hydrology}{\gdef\@journal{hydrology}}
\DeclareOption{informatics}{\gdef\@journal{informatics}}
\DeclareOption{information}{\gdef\@journal{information}}
\DeclareOption{infrastructures}{\gdef\@journal{infrastructures}}
\DeclareOption{inorganics}{\gdef\@journal{inorganics}}
\DeclareOption{insects}{\gdef\@journal{insects}}
\DeclareOption{instruments}{\gdef\@journal{instruments}}
\DeclareOption{ijerph}{\gdef\@journal{ijerph}}
\DeclareOption{ijfs}{\gdef\@journal{ijfs}}
\DeclareOption{ijms}{\gdef\@journal{ijms}}
\DeclareOption{ijgi}{\gdef\@journal{ijgi}}
\DeclareOption{inventions}{\gdef\@journal{inventions}}
\DeclareOption{jcdd}{\gdef\@journal{jcdd}}
\DeclareOption{jcm}{\gdef\@journal{jcm}}
\DeclareOption{jdb}{\gdef\@journal{jdb}}
\DeclareOption{jfb}{\gdef\@journal{jfb}}
\DeclareOption{jfmk}{\gdef\@journal{jfmk}}
\DeclareOption{jof}{\gdef\@journal{jof}}
\DeclareOption{jimaging}{\gdef\@journal{jimaging}}
\DeclareOption{jintelligence}{\gdef\@journal{jintelligence}}
\DeclareOption{jlpea}{\gdef\@journal{jlpea}}
\DeclareOption{jmse}{\gdef\@journal{jmse}}
\DeclareOption{jpm}{\gdef\@journal{jpm}}
\DeclareOption{jrfm}{\gdef\@journal{jrfm}}
\DeclareOption{jsan}{\gdef\@journal{jsan}}
\DeclareOption{land}{\gdef\@journal{land}}
\DeclareOption{languages}{\gdef\@journal{languages}}
\DeclareOption{laws}{\gdef\@journal{laws}}
\DeclareOption{life}{\gdef\@journal{life}}
\DeclareOption{literature}{\gdef\@journal{literature}}
\DeclareOption{lubricants}{\gdef\@journal{lubricants}}
\DeclareOption{machines}{\gdef\@journal{machines}}
\DeclareOption{magnetochemistry}{\gdef\@journal{magnetochemistry}}
\DeclareOption{marinedrugs}{\gdef\@journal{marinedrugs}}
\DeclareOption{materials}{\gdef\@journal{materials}}
\DeclareOption{mathematics}{\gdef\@journal{mathematics}}
\DeclareOption{mca}{\gdef\@journal{mca}}
\DeclareOption{medsci}{\gdef\@journal{medsci}}
\DeclareOption{medicalsciences}{\gdef\@journal{medsci}} %to have backwards compatibility 
\DeclareOption{medicines}{\gdef\@journal{medicines}}
\DeclareOption{membranes}{\gdef\@journal{membranes}}
\DeclareOption{metabolites}{\gdef\@journal{metabolites}}
\DeclareOption{metals}{\gdef\@journal{metals}}
\DeclareOption{microarrays}{\gdef\@journal{microarrays}}
\DeclareOption{micromachines}{\gdef\@journal{micromachines}}
\DeclareOption{microorganisms}{\gdef\@journal{microorganisms}}
\DeclareOption{minerals}{\gdef\@journal{minerals}}
\DeclareOption{molbank}{\gdef\@journal{molbank}}
\DeclareOption{molecules}{\gdef\@journal{molecules}}
\DeclareOption{mps}{\gdef\@journal{mps}}
\DeclareOption{nanomaterials}{\gdef\@journal{nanomaterials}}
\DeclareOption{ncrna}{\gdef\@journal{ncrna}}
\DeclareOption{neonatalscreening}{\gdef\@journal{neonatalscreening}}
\DeclareOption{nutrients}{\gdef\@journal{nutrients}}
\DeclareOption{particles}{\gdef\@journal{particles}}
\DeclareOption{pathogens}{\gdef\@journal{pathogens}}
\DeclareOption{pharmaceuticals}{\gdef\@journal{pharmaceuticals}}
\DeclareOption{pharmaceutics}{\gdef\@journal{pharmaceutics}}
\DeclareOption{pharmacy}{\gdef\@journal{pharmacy}}
\DeclareOption{philosophies}{\gdef\@journal{philosophies}}
\DeclareOption{photonics}{\gdef\@journal{photonics}}
\DeclareOption{plants}{\gdef\@journal{plants}}
\DeclareOption{polymers}{\gdef\@journal{polymers}}
\DeclareOption{processes}{\gdef\@journal{processes}}
\DeclareOption{proteomes}{\gdef\@journal{proteomes}}
\DeclareOption{publications}{\gdef\@journal{publications}}
\DeclareOption{recycling}{\gdef\@journal{recycling}}
\DeclareOption{religions}{\gdef\@journal{religions}}
\DeclareOption{remotesensing}{\gdef\@journal{remotesensing}}
\DeclareOption{resources}{\gdef\@journal{resources}}
\DeclareOption{risks}{\gdef\@journal{risks}}
\DeclareOption{robotics}{\gdef\@journal{robotics}}
\DeclareOption{safety}{\gdef\@journal{safety}}
\DeclareOption{sensors}{\gdef\@journal{sensors}}
\DeclareOption{sexes}{\gdef\@journal{sexes}}
\DeclareOption{sinusitis}{\gdef\@journal{sinusitis}}
\DeclareOption{socsci}{\gdef\@journal{socsci}}
\DeclareOption{socialsciences}{\gdef\@journal{socsci}} %to have backwards compatibility 
\DeclareOption{societies}{\gdef\@journal{societies}}
\DeclareOption{soils}{\gdef\@journal{soils}}
\DeclareOption{sports}{\gdef\@journal{sports}}
\DeclareOption{standards}{\gdef\@journal{standards}}
\DeclareOption{sustainability}{\gdef\@journal{sustainability}}
\DeclareOption{symmetry}{\gdef\@journal{symmetry}}
\DeclareOption{systems}{\gdef\@journal{systems}}
\DeclareOption{technologies}{\gdef\@journal{technologies}}
\DeclareOption{toxics}{\gdef\@journal{toxics}}
\DeclareOption{toxins}{\gdef\@journal{toxins}}
\DeclareOption{universe}{\gdef\@journal{universe}}
\DeclareOption{urbansci}{\gdef\@journal{urbansci}}
\DeclareOption{vaccines}{\gdef\@journal{vaccines}}
\DeclareOption{vetsci}{\gdef\@journal{vetsci}}
\DeclareOption{veterinarysciences}{\gdef\@journal{vetsci}} %to have backwards compatibility
\DeclareOption{viruses}{\gdef\@journal{viruses}}
\DeclareOption{visions}{\gdef\@journal{visions}}
\DeclareOption{water}{\gdef\@journal{water}}

%% Uppercase will have the same effect
\DeclareOption{Actuators}{\gdef\@journal{actuators}}
\DeclareOption{Admsci}{\gdef\@journal{admsci}}
\DeclareOption{Administrativesciences}{\gdef\@journal{admsci}} %to have backwards compatibility
\DeclareOption{Aerospace}{\gdef\@journal{aerospace}}
\DeclareOption{Agriculture}{\gdef\@journal{agriculture}}
\DeclareOption{Agronomy}{\gdef\@journal{agronomy}}
\DeclareOption{Algorithms}{\gdef\@journal{algorithms}}
\DeclareOption{Animals}{\gdef\@journal{animals}}
\DeclareOption{Antibiotics}{\gdef\@journal{antibiotics}}
\DeclareOption{Antibodies}{\gdef\@journal{antibodies}}
\DeclareOption{Antioxidants}{\gdef\@journal{antioxidants}}
\DeclareOption{Applsci}{\gdef\@journal{applsci}}
\DeclareOption{Appliedsciences}{\gdef\@journal{applsci}} %to have backwards compatibility
\DeclareOption{Arts}{\gdef\@journal{arts}}
\DeclareOption{Atmosphere}{\gdef\@journal{atmosphere}}
\DeclareOption{Atoms}{\gdef\@journal{atoms}}
\DeclareOption{Axioms}{\gdef\@journal{axioms}}
\DeclareOption{Batteries}{\gdef\@journal{batteries}}
\DeclareOption{Behavsci}{\gdef\@journal{behavsci}}
\DeclareOption{Behavioralsciences}{\gdef\@journal{behavsci}} %to have backwards compatibility
\DeclareOption{Beverages}{\gdef\@journal{beverages}}
\DeclareOption{Bioengineering}{\gdef\@journal{bioengineering}}
\DeclareOption{Biology}{\gdef\@journal{biology}}
\DeclareOption{Biomedicines}{\gdef\@journal{biomedicines}}
\DeclareOption{Biomimetics}{\gdef\@journal{biomimetics}}
\DeclareOption{Biomolecules}{\gdef\@journal{biomolecules}}
\DeclareOption{Biosensors}{\gdef\@journal{biosensors}}
\DeclareOption{Brainsci}{\gdef\@journal{brainsci}}
\DeclareOption{Brainsciences}{\gdef\@journal{brainsci}} %to have backwards compatibility
\DeclareOption{Buildings}{\gdef\@journal{buildings}}
\DeclareOption{C}{\gdef\@journal{carbon}}
\DeclareOption{Carbon}{\gdef\@journal{carbon}}
\DeclareOption{Cancers}{\gdef\@journal{cancers}}
\DeclareOption{Catalysts}{\gdef\@journal{catalysts}}
\DeclareOption{Cells}{\gdef\@journal{cells}}
\DeclareOption{Challenges}{\gdef\@journal{challenges}}
\DeclareOption{Chemosensors}{\gdef\@journal{chemosensors}}
\DeclareOption{Children}{\gdef\@journal{children}}
\DeclareOption{Chromatography}{\gdef\@journal{chromatography}}
\DeclareOption{Climate}{\gdef\@journal{climate}}
\DeclareOption{Coatings}{\gdef\@journal{coatings}}
\DeclareOption{Computation}{\gdef\@journal{computation}}
\DeclareOption{Computers}{\gdef\@journal{computers}}
\DeclareOption{Condensedmatter}{\gdef\@journal{condensedmatter}}
\DeclareOption{Cosmetics}{\gdef\@journal{cosmetics}}
\DeclareOption{Cryptography}{\gdef\@journal{cryptography}}
\DeclareOption{Crystals}{\gdef\@journal{crystals}}
\DeclareOption{Data}{\gdef\@journal{data}}
\DeclareOption{Dentistry}{\gdef\@journal{dentistry}}
\DeclareOption{Dentistryjournal}{\gdef\@journal{dentistry}} %to have backwards compatibility 
\DeclareOption{Designs}{\gdef\@journal{designs}}
\DeclareOption{Diagnostics}{\gdef\@journal{diagnostics}}
\DeclareOption{Diseases}{\gdef\@journal{diseases}}
\DeclareOption{Diversity}{\gdef\@journal{diversity}}
\DeclareOption{Econometrics}{\gdef\@journal{econometrics}}
\DeclareOption{Economies}{\gdef\@journal{economies}}
\DeclareOption{Education}{\gdef\@journal{education}}
\DeclareOption{Educationsciences}{\gdef\@journal{education}} %to have backwards compatibility
\DeclareOption{Electronics}{\gdef\@journal{electronics}}
\DeclareOption{Energies}{\gdef\@journal{energies}}
\DeclareOption{Entropy}{\gdef\@journal{entropy}}
\DeclareOption{Environments}{\gdef\@journal{environments}}
\DeclareOption{Epigenomes}{\gdef\@journal{epigenomes}}
\DeclareOption{Fermentation}{\gdef\@journal{fermentation}}
\DeclareOption{Fibers}{\gdef\@journal{fibers}}
\DeclareOption{Fishes}{\gdef\@journal{fishes}}
\DeclareOption{Fluids}{\gdef\@journal{fluids}}
\DeclareOption{Foods}{\gdef\@journal{foods}}
\DeclareOption{Forests}{\gdef\@journal{forests}}
\DeclareOption{Futureinternet}{\gdef\@journal{futureinternet}}
\DeclareOption{Galaxies}{\gdef\@journal{galaxies}}
\DeclareOption{Games}{\gdef\@journal{games}}
\DeclareOption{Gels}{\gdef\@journal{gels}}
\DeclareOption{Genealogy}{\gdef\@journal{genealogy}}
\DeclareOption{Genes}{\gdef\@journal{genes}}
\DeclareOption{Geosciences}{\gdef\@journal{geosciences}}
\DeclareOption{Geriatrics}{\gdef\@journal{geriatrics}}
\DeclareOption{Healthcare}{\gdef\@journal{healthcare}}
\DeclareOption{Horticulturae}{\gdef\@journal{horticulturae}}
\DeclareOption{Humanities}{\gdef\@journal{humanities}}
\DeclareOption{Hydrology}{\gdef\@journal{hydrology}}
\DeclareOption{Informatics}{\gdef\@journal{informatics}}
\DeclareOption{Information}{\gdef\@journal{information}}
\DeclareOption{Infrastructures}{\gdef\@journal{infrastructures}}
\DeclareOption{Inorganics}{\gdef\@journal{inorganics}}
\DeclareOption{Insects}{\gdef\@journal{insects}}
\DeclareOption{Instruments}{\gdef\@journal{instruments}}
\DeclareOption{Ijerph}{\gdef\@journal{ijerph}}
\DeclareOption{Ijfs}{\gdef\@journal{ijfs}}
\DeclareOption{Ijms}{\gdef\@journal{ijms}}
\DeclareOption{Ijgi}{\gdef\@journal{ijgi}}
\DeclareOption{Inventions}{\gdef\@journal{inventions}}
\DeclareOption{Jcdd}{\gdef\@journal{jcdd}}
\DeclareOption{Jcm}{\gdef\@journal{jcm}}
\DeclareOption{Jdb}{\gdef\@journal{jdb}}
\DeclareOption{Jfb}{\gdef\@journal{jfb}}
\DeclareOption{Jmk}{\gdef\@journal{jfmk}}
\DeclareOption{Jof}{\gdef\@journal{jof}}
\DeclareOption{Jimaging}{\gdef\@journal{jimaging}}
\DeclareOption{Jintelligence}{\gdef\@journal{jintelligence}}
\DeclareOption{Jlpea}{\gdef\@journal{jlpea}}
\DeclareOption{Jmse}{\gdef\@journal{jmse}}
\DeclareOption{Jpm}{\gdef\@journal{jpm}}
\DeclareOption{Jrfm}{\gdef\@journal{jrfm}}
\DeclareOption{Jsan}{\gdef\@journal{jsan}}
\DeclareOption{Land}{\gdef\@journal{land}}
\DeclareOption{Languages}{\gdef\@journal{languages}}
\DeclareOption{Laws}{\gdef\@journal{laws}}
\DeclareOption{Life}{\gdef\@journal{life}}
\DeclareOption{Literature}{\gdef\@journal{literature}}
\DeclareOption{Lubricants}{\gdef\@journal{lubricants}}
\DeclareOption{Machines}{\gdef\@journal{machines}}
\DeclareOption{Magnetochemistry}{\gdef\@journal{magnetochemistry}}
\DeclareOption{MarineDrugs}{\gdef\@journal{marinedrugs}}
\DeclareOption{Materials}{\gdef\@journal{materials}}
\DeclareOption{Mathematics}{\gdef\@journal{mathematics}}
\DeclareOption{Mca}{\gdef\@journal{mca}}
\DeclareOption{Medsci}{\gdef\@journal{medsci}}
\DeclareOption{Medicalsciences}{\gdef\@journal{medsci}} %to have backwards compatibility
\DeclareOption{Medicines}{\gdef\@journal{medicines}}
\DeclareOption{Membranes}{\gdef\@journal{membranes}}
\DeclareOption{Metabolites}{\gdef\@journal{metabolites}}
\DeclareOption{Metals}{\gdef\@journal{metals}}
\DeclareOption{Microarrays}{\gdef\@journal{microarrays}}
\DeclareOption{Micromachines}{\gdef\@journal{micromachines}}
\DeclareOption{Microorganisms}{\gdef\@journal{microorganisms}}
\DeclareOption{Minerals}{\gdef\@journal{minerals}}
\DeclareOption{Molbank}{\gdef\@journal{molbank}}
\DeclareOption{Molecules}{\gdef\@journal{molecules}}
\DeclareOption{Mps}{\gdef\@journal{mps}}
\DeclareOption{Nanomaterials}{\gdef\@journal{nanomaterials}}
\DeclareOption{ncRNA}{\gdef\@journal{ncrna}}
\DeclareOption{Neonatalscreening}{\gdef\@journal{neonatalscreening}}
\DeclareOption{Nutrients}{\gdef\@journal{nutrients}}
\DeclareOption{Particles}{\gdef\@journal{particles}}
\DeclareOption{Pathogens}{\gdef\@journal{pathogens}}
\DeclareOption{Pharmaceuticals}{\gdef\@journal{pharmaceuticals}}
\DeclareOption{Pharmaceutics}{\gdef\@journal{pharmaceutics}}
\DeclareOption{Pharmacy}{\gdef\@journal{pharmacy}}
\DeclareOption{Philosophies}{\gdef\@journal{philosophies}}
\DeclareOption{Photonics}{\gdef\@journal{photonics}}
\DeclareOption{Plants}{\gdef\@journal{plants}}
\DeclareOption{Polymers}{\gdef\@journal{polymers}}
\DeclareOption{Processes}{\gdef\@journal{processes}}
\DeclareOption{Proteomes}{\gdef\@journal{proteomes}}
\DeclareOption{Publications}{\gdef\@journal{publications}}
\DeclareOption{Recycling}{\gdef\@journal{recycling}}
\DeclareOption{Religions}{\gdef\@journal{religions}}
\DeclareOption{RemoteSensing}{\gdef\@journal{remotesensing}}
\DeclareOption{Resources}{\gdef\@journal{resources}}
\DeclareOption{Risks}{\gdef\@journal{risks}}
\DeclareOption{Robotics}{\gdef\@journal{robotics}}
\DeclareOption{Safety}{\gdef\@journal{safety}}
\DeclareOption{Sensors}{\gdef\@journal{sensors}}
\DeclareOption{Sexes}{\gdef\@journal{sexes}}
\DeclareOption{Sinusitis}{\gdef\@journal{sinusitis}}
\DeclareOption{Socsci}{\gdef\@journal{socsci}}
\DeclareOption{SocialSciences}{\gdef\@journal{socsci}} %to have backwards compatibility
\DeclareOption{Societies}{\gdef\@journal{societies}}
\DeclareOption{Soils}{\gdef\@journal{soils}}
\DeclareOption{Sports}{\gdef\@journal{sports}}
\DeclareOption{Standards}{\gdef\@journal{standards}}
\DeclareOption{Sustainability}{\gdef\@journal{sustainability}}
\DeclareOption{Symmetry}{\gdef\@journal{symmetry}}
\DeclareOption{Systems}{\gdef\@journal{systems}}
\DeclareOption{Technologies}{\gdef\@journal{technologies}}
\DeclareOption{Toxics}{\gdef\@journal{toxics}}
\DeclareOption{Toxins}{\gdef\@journal{toxins}}
\DeclareOption{Universe}{\gdef\@journal{universe}}
\DeclareOption{Urbansci}{\gdef\@journal{urbansci}}
\DeclareOption{Vaccines}{\gdef\@journal{vaccines}}
\DeclareOption{Vetsci}{\gdef\@journal{vetsci}}
\DeclareOption{VeterinarySciences}{\gdef\@journal{vetsci}} %to have backwards compatibility
\DeclareOption{Viruses}{\gdef\@journal{viruses}}
\DeclareOption{Visions}{\gdef\@journal{visions}}
\DeclareOption{Water}{\gdef\@journal{water}}
\DeclareOption{journal}{\ClassWarning{mdpi}{You have not specified the journal. The first class option of the document class mdpi specifies the journal. The word journal should be replaced by one of the following journal names: actuators, admsci, aerospace, agriculture, agronomy, algorithms, animals, antibiotics, antibodies, antioxidants, applsci, arts, atmosphere, atoms, axioms, batteries, behavsci, beverages, bioengineering, biology, biomedicines, biomimetics, biomolecules, biosensors, brainsci, buildings, carbon, cancers, catalysts, cells, challenges, chemosensors, children, chromatography, climate, coatings, computation, computers, condensedmatter, cosmetics, cryptography, crystals, data, dentistry, designs, diagnostics, diseases, diversity, econometrics, economies, education, electronics, energies, entropy, environments, epigenomes, fermentation, fibers, fishes, fluids, foods, forests, futureinternet, galaxies, games, gels, genealogy, genes, geosciences, geriatrics, healthcare, horticulturae, humanities, hydrology, informatics, information, infrastructures, inorganics, insects, instruments, ijerph, ijfs, ijms, ijgi, inventions, jcdd, jcm, jdb, jfb, jfmk, jimaging, jof, jintelligence, jlpea, jmse, jpm, jrfm, jsan, land, languages, laws, life, literature, lubricants, machines, magnetochemistry, marinedrugs, materials, mathematics, mca, medsci, medicines, membranes, metabolites, metals, microarrays, micromachines, microorganisms, minerals, molbank, molecules, mps, nanomaterials, ncrna, neonatalscreening, nutrients, particles, pathogens, pharmaceuticals, pharmaceutics, pharmacy, philosophies, photonics, plants, polymers, processes, proteomes, publications, recycling, religions, remotesensing, resources, risks, robotics, safety, sensors, sexes, sinusitis, socsci, societies, soils, sports, standards, sustainability, symmetry, systems, technologies, toxics, toxins, universe, urbansci, vaccines, vetsci, viruses, visions, water}} 

%% To choose the type of manuscript
\DeclareOption{addendum}{\gdef\@arttype{Addendum}}
\DeclareOption{article}{\gdef\@arttype{Article}}
\DeclareOption{bookreview}{\gdef\@arttype{Book Review}}
\DeclareOption{briefreport}{\gdef\@arttype{Brief Report}}
\DeclareOption{casereport}{\gdef\@arttype{Case Report}}
\DeclareOption{changes}{\gdef\@arttype{Changes}}
\DeclareOption{comment}{\gdef\@arttype{Comment}}
\DeclareOption{commentary}{\gdef\@arttype{Commentary}}
\DeclareOption{communication}{\gdef\@arttype{Communication}}
\DeclareOption{conceptpaper}{\gdef\@arttype{Concept Paper}}
\DeclareOption{correction}{\gdef\@arttype{Correction}}
\DeclareOption{conferencereport}{\gdef\@arttype{Conference Report}}
\DeclareOption{meetingreport}{\gdef\@arttype{Meeting Report}}
\DeclareOption{creative}{\gdef\@arttype{Creative}}
\DeclareOption{datadescriptor}{\gdef\@arttype{Data Descriptor}}
\DeclareOption{discussion}{\gdef\@arttype{Discussion}}
\DeclareOption{editorial}{\gdef\@arttype{Editorial}}
\DeclareOption{essay}{\gdef\@arttype{Essay}}
\DeclareOption{erratum}{\gdef\@arttype{Erratum}}
\DeclareOption{hypothesis}{\gdef\@arttype{Hypothesis}}
\DeclareOption{interestingimage}{\gdef\@arttype{Interesting Image}}
\DeclareOption{letter}{\gdef\@arttype{Letter}}
\DeclareOption{newbookreceived}{\gdef\@arttype{New Book Received}}
\DeclareOption{obituary}{\gdef\@arttype{Obituary}}
\DeclareOption{opinion}{\gdef\@arttype{Opinion}}
\DeclareOption{projectreport}{\gdef\@arttype{Project Report}}
\DeclareOption{reply}{\gdef\@arttype{Reply}}
\DeclareOption{retraction}{\gdef\@arttype{Retraction}}
\DeclareOption{review}{\gdef\@arttype{Review}}
\DeclareOption{shortnote}{\gdef\@arttype{Short Note}}
\DeclareOption{supfile}{\gdef\@arttype{Supfile}}
\DeclareOption{technicalnote}{\gdef\@arttype{Technical Note}}

%% To choose the status of the manuscript
\DeclareOption{submit}{\gdef\@status{submit}}
\DeclareOption{accept}{\gdef\@status{accept}}

%% To choose the whether there is one or more authors
\DeclareOption{oneauthor}{\gdef\@authornum{author}}
\DeclareOption{moreauthors}{\gdef\@authornum{authors}}

%% Add the chosen options to the class
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}

%% Defaults
\ExecuteOptions{entropy,10pt,a4paper,article,submit,oneauthor}

%% Process options
\ProcessOptions\relax

%% MORE DECLARATIONS

%%%% Maths environments
\RequirePackage{amsthm}
\newtheoremstyle{mdpi}% name
{12pt}% space above
{12pt}% space below
{\itshape}% body font
{}% indent amount 1
{\bfseries}% theorem head font
{.}% punctuation after theorem head
{.5em}% space after theorem head
{}% theorem head spec (can be left empty, meaning `normal')

\newtheoremstyle{mdpidefinition}% name
{12pt}% space above
{12pt}% space below
{}% body font
{}% hndent amount 1
{\bfseries}% theorem head font
{.}% punctuation after theorem head
{.5em}% space after theorem head
{}% theorem head spec (can be left empty, meaning `normal')

\renewcommand{\qed}{\unskip\nobreak\quad\qedsymbol} %% This places the symbol right after the text instead of placing it at the end on the line.

\renewenvironment{proof}[1][\proofname]{\par %% \proofname allows to have "Proof of my theorem"
  \pushQED{\qed}%
  \normalfont \topsep6\p@\@plus6\p@\relax
  \trivlist
  \item[\hskip\labelsep
        \bfseries %% "Proof" is bold
    #1\@addpunct{.}]\ignorespaces %% Period instead of colon
}{%
  \popQED\endtrivlist\@endpefalse
}

%%%% Hyphenation
\RequirePackage[none]{hyphenat}
\sloppy

%%%% References
\RequirePackage[sort&compress]{natbib}
\bibpunct{[}{]}{,}{n}{}{}

\renewcommand\NAT@set@cites{%
  \ifNAT@numbers
    \ifNAT@super \let\@cite\NAT@citesuper
       \def\NAT@mbox##1{\unskip\nobreak\textsuperscript{##1}}%
       \let\citeyearpar=\citeyear
       \let\NAT@space\relax
       \def\NAT@super@kern{\kern\p@}%
    \else
       \let\NAT@mbox=\mbox
       \let\@cite\NAT@citenum
       \let\NAT@space\relax
       \let\NAT@super@kern\relax
    \fi
    \let\@citex\NAT@citexnum
    \let\@biblabel\NAT@biblabelnum
    \let\@bibsetup\NAT@bibsetnum
    \renewcommand\NAT@idxtxt{\NAT@name\NAT@spacechar\NAT@open\NAT@num\NAT@close}%
    \def\natexlab##1{}%
    \def\NAT@penalty{\penalty\@m}%
  \else
    \let\@cite\NAT@cite
    \let\@citex\NAT@citex
    \let\@biblabel\NAT@biblabel
    \let\@bibsetup\NAT@bibsetup
    \let\NAT@space\NAT@spacechar
    \let\NAT@penalty\@empty
    \renewcommand\NAT@idxtxt{\NAT@name\NAT@spacechar\NAT@open\NAT@date\NAT@close}%
    \def\natexlab##1{##1}%
  \fi}

%%%%% Hyperlinks
\RequirePackage[unicode=true,
bookmarksopen={true},
pdffitwindow=true, 
colorlinks=true, 
linkcolor=blue, 
citecolor=blue, 
urlcolor=blue, 
hyperfootnotes=false, 
pdfstartview={FitH},
pdfpagemode= None]{hyperref}   

%% To have the option to change the urlcolor
\newcommand{\changeurlcolor}[1]{\hypersetup{urlcolor=#1}} 

%% Metadata
\newcommand{\org@maketitle}{}% LATEX-Check
\let\org@maketitle\maketitle
\def\maketitle{%
\hypersetup{
pdftitle={\@Title},
pdfsubject={\@abstract},
pdfkeywords={\@keyword},
pdfauthor={\@AuthorNames}
}%
\org@maketitle
}

%%%% Footnotes
\RequirePackage[hang]{footmisc}
\setlength{\skip\footins}{1.2cm}
\setlength{\footnotemargin}{5mm}
\def\footnoterule{\kern-14\p@
\hrule \@width 2in \kern 11.6\p@}

%%%% URL
\RequirePackage{url}
\urlstyle{same}

%%%% Widows & orphans
\clubpenalty=10000
\widowpenalty=10000
\displaywidowpenalty=10000

%%%% Front matter
\newcommand{\firstargument}{}
\newcommand{\Title}[1]{\gdef\@Title{#1}}%
\newcommand{\Author}[1]{\gdef\@Author{#1}}%
\def\@AuthorNames{}
\newcommand{\AuthorNames}[1]{\gdef\@AuthorNames{#1}}%
\newcommand{\firstpage}[1]{\gdef\@firstpage{#1}}
\newcommand{\doinum}[1]{\gdef\@doinum{#1}}
\newcommand{\pubvolume}[1]{\gdef\@pubvolume{#1}}
\newcommand{\pubyear}[1]{\gdef\@pubyear{#1}}
\newcommand{\copyrightyear}[1]{\gdef\@copyrightyear{#1}}
\newcommand{\address}[2][]{\renewcommand{\firstargument}{#1}\gdef\@address{#2}}
\newcommand{\corresfirstargument}{}
\def\@corres{}
\newcommand{\corres}[2][]{\renewcommand{\corresfirstargument}{#1}\gdef\@corres{#2}} 
\def\@conference{}
\newcommand{\conference}[1]{\gdef\@conference{#1}}%
\def\@abstract{}
\renewcommand{\abstract}[1]{\gdef\@abstract{#1}}
\def\@externaleditor{}
\newcommand{\externaleditor}[1]{\gdef\@externaleditor{#1}}
\newcommand{\history}[1]{\gdef\@history{#1}} 
\def\@pacs{}
\newcommand{\PACS}[1]{\gdef\@pacs{#1}} 
\def\@msc{}
\newcommand{\MSC}[1]{\gdef\@msc{#1}} 
\def\@jel{}
\newcommand{\JEL}[1]{\gdef\@jel{#1}}
\def\@keyword{}
\newcommand{\keyword}[1]{\gdef\@keyword{#1}}
\def\@dataset{}
\newcommand{\dataset}[1]{\gdef\@dataset{#1}}
\def\@datasetlicense{}
\newcommand{\datasetlicense}[1]{\gdef\@datasetlicense{#1}}
\def\@firstnote{}
\newcommand{\firstnote}[1]{\gdef\@firstnote{#1}}
\def\@secondnote{}
\newcommand{\secondnote}[1]{\gdef\@secondnote{#1}}%
\def\@simplesumm{}
\newcommand{\simplesumm}[1]{\gdef\@simplesumm{#1}}
\newcommand{\articlenumber}[1]{\gdef\@articlenumber{#1}}

%%%% Journal name for the header
\newcommand{\journalname}{%
  \ifthenelse{\equal{\@journal}{actuators}}{Actuators}{}%
  \ifthenelse{\equal{\@journal}{admsci}}{Adm. Sci.}{}%
  \ifthenelse{\equal{\@journal}{aerospace}}{Aerospace}{}%
  \ifthenelse{\equal{\@journal}{agriculture}}{Agriculture}{}%
  \ifthenelse{\equal{\@journal}{agronomy}}{Agronomy}{}%
  \ifthenelse{\equal{\@journal}{algorithms}}{Algorithms}{}%
  \ifthenelse{\equal{\@journal}{animals}}{Animals}{}%
  \ifthenelse{\equal{\@journal}{antibiotics}}{Antibiotics}{}%
  \ifthenelse{\equal{\@journal}{antibodies}}{Antibodies}{}%
  \ifthenelse{\equal{\@journal}{antioxidants}}{Antioxidants}{}%
  \ifthenelse{\equal{\@journal}{applsci}}{Appl. Sci.}{}%
  \ifthenelse{\equal{\@journal}{arts}}{Arts}{}%
  \ifthenelse{\equal{\@journal}{atmosphere}}{Atmosphere}{}%
  \ifthenelse{\equal{\@journal}{atoms}}{Atoms}{}%
  \ifthenelse{\equal{\@journal}{axioms}}{Axioms}{}%
  \ifthenelse{\equal{\@journal}{batteries}}{Batteries}{}%
  \ifthenelse{\equal{\@journal}{behavsci}}{Behav. Sci.}{}%
  \ifthenelse{\equal{\@journal}{beverages}}{Beverages}{}%
  \ifthenelse{\equal{\@journal}{bioengineering}}{Bioengineering}{}%
  \ifthenelse{\equal{\@journal}{biology}}{Biology}{}%
  \ifthenelse{\equal{\@journal}{biomedicines}}{Biomedicines}{}%
  \ifthenelse{\equal{\@journal}{biomimetics}}{Biomimetics}{}%
  \ifthenelse{\equal{\@journal}{biomolecules}}{Biomolecules}{}%
  \ifthenelse{\equal{\@journal}{biosensors}}{Biosensors}{}%
  \ifthenelse{\equal{\@journal}{brainsci}}{Brain Sci.}{}%
  \ifthenelse{\equal{\@journal}{buildings}}{Buildings}{}%
  \ifthenelse{\equal{\@journal}{carbon}}{C}{}%
  \ifthenelse{\equal{\@journal}{cancers}}{Cancers}{}%
  \ifthenelse{\equal{\@journal}{catalysts}}{Catalysts}{}%
  \ifthenelse{\equal{\@journal}{cells}}{Cells}{}%
  \ifthenelse{\equal{\@journal}{challenges}}{Challenges}{}%
  \ifthenelse{\equal{\@journal}{chemosensors}}{Chemosensors}{}%
  \ifthenelse{\equal{\@journal}{children}}{Children}{}%
  \ifthenelse{\equal{\@journal}{chromatography}}{Chromatography}{}%
  \ifthenelse{\equal{\@journal}{climate}}{Climate}{}%
  \ifthenelse{\equal{\@journal}{coatings}}{Coatings}{}%
  \ifthenelse{\equal{\@journal}{computation}}{Computation}{}%
  \ifthenelse{\equal{\@journal}{computers}}{Computers}{}%
  \ifthenelse{\equal{\@journal}{condensedmatter}}{Condens. Matter}{}%
  \ifthenelse{\equal{\@journal}{cosmetics}}{Cosmetics}{}%
  \ifthenelse{\equal{\@journal}{cryptography}}{Cryptography}{}%
  \ifthenelse{\equal{\@journal}{crystals}}{Crystals}{}%
  \ifthenelse{\equal{\@journal}{data}}{Data}{}%
  \ifthenelse{\equal{\@journal}{dentistry}}{Dent. J.}{}%
  \ifthenelse{\equal{\@journal}{designs}}{Designs}{}%
  \ifthenelse{\equal{\@journal}{diagnostics}}{Diagnostics}{}%
  \ifthenelse{\equal{\@journal}{diseases}}{Diseases}{}%
  \ifthenelse{\equal{\@journal}{diversity}}{Diversity}{}%
  \ifthenelse{\equal{\@journal}{econometrics}}{Econometrics}{}%
  \ifthenelse{\equal{\@journal}{economies}}{Economies}{}%
  \ifthenelse{\equal{\@journal}{education}}{Educ. Sci.}{}%
  \ifthenelse{\equal{\@journal}{electronics}}{Electronics}{}%
  \ifthenelse{\equal{\@journal}{energies}}{Energies}{}%
  \ifthenelse{\equal{\@journal}{entropy}}{Entropy}{}%
  \ifthenelse{\equal{\@journal}{environments}}{Environments}{}%
  \ifthenelse{\equal{\@journal}{epigenomes}}{Epigenomes}{}%
  \ifthenelse{\equal{\@journal}{fermentation}}{Fermentation}{}%
  \ifthenelse{\equal{\@journal}{fibers}}{Fibers}{}%
  \ifthenelse{\equal{\@journal}{fishes}}{Fishes}{}%
  \ifthenelse{\equal{\@journal}{fluids}}{Fluids}{}%
  \ifthenelse{\equal{\@journal}{foods}}{Foods}{}%
  \ifthenelse{\equal{\@journal}{forests}}{Forests}{}%
  \ifthenelse{\equal{\@journal}{futureinternet}}{Future Internet}{}%
  \ifthenelse{\equal{\@journal}{galaxies}}{Galaxies}{}%
  \ifthenelse{\equal{\@journal}{games}}{Games}{}%
  \ifthenelse{\equal{\@journal}{gels}}{Gels}{}%
  \ifthenelse{\equal{\@journal}{genealogy}}{Genealogy}{}%
  \ifthenelse{\equal{\@journal}{genes}}{Genes}{}%
  \ifthenelse{\equal{\@journal}{geosciences}}{Geosciences}{}%
  \ifthenelse{\equal{\@journal}{geriatrics}}{Geriatrics}{}%
  \ifthenelse{\equal{\@journal}{healthcare}}{Healthcare}{}%
  \ifthenelse{\equal{\@journal}{horticulturae}}{Horticulturae}{}%
  \ifthenelse{\equal{\@journal}{humanities}}{Humanities}{}%
  \ifthenelse{\equal{\@journal}{hydrology}}{Hydrology}{}%
  \ifthenelse{\equal{\@journal}{informatics}}{Informatics}{}%
  \ifthenelse{\equal{\@journal}{information}}{Information}{}%
  \ifthenelse{\equal{\@journal}{infrastructures}}{Infrastructures}{}%
  \ifthenelse{\equal{\@journal}{inorganics}}{Inorganics}{}%
  \ifthenelse{\equal{\@journal}{insects}}{Insects}{}%
  \ifthenelse{\equal{\@journal}{instruments}}{Instruments}{}%
  \ifthenelse{\equal{\@journal}{ijerph}}{Int. J. Environ. Res. Public Health}{}%
  \ifthenelse{\equal{\@journal}{ijfs}}{Int. J. Financial Stud.}{}%
  \ifthenelse{\equal{\@journal}{ijms}}{Int. J. Mol. Sci.}{}%
  \ifthenelse{\equal{\@journal}{ijgi}}{ISPRS Int. J. Geo-Inf.}{}%
  \ifthenelse{\equal{\@journal}{inventions}}{Inventions}{}%
  \ifthenelse{\equal{\@journal}{jcdd}}{J. Cardiovasc. Dev. Dis.}{}%
  \ifthenelse{\equal{\@journal}{jcm}}{J. Clin. Med.}{}%
  \ifthenelse{\equal{\@journal}{jdb}}{J. Dev. Biol.}{}%
  \ifthenelse{\equal{\@journal}{jfb}}{J. Funct. Biomater.}{}%   
  \ifthenelse{\equal{\@journal}{jfmk}}{J. Funct. Morphol. Kinesiol.}{}%   
  \ifthenelse{\equal{\@journal}{jof}}{J. Fungi}{}%
  \ifthenelse{\equal{\@journal}{jimaging}}{J. Imaging}{}%
  \ifthenelse{\equal{\@journal}{jintelligence}}{J. Intell.}{}% 
  \ifthenelse{\equal{\@journal}{jlpea}}{J. Low Power Electron. Appl.}{}%
  \ifthenelse{\equal{\@journal}{jmse}}{J. Mar. Sci. Eng.}{}%??????
  \ifthenelse{\equal{\@journal}{jpm}}{J. Pers. Med.}{}%
  \ifthenelse{\equal{\@journal}{jrfm}}{J. Risk Financial Manag.}{}%
  \ifthenelse{\equal{\@journal}{jsan}}{J. Sens. Actuator Netw.}{}%
  \ifthenelse{\equal{\@journal}{land}}{Land}{}%
  \ifthenelse{\equal{\@journal}{languages}}{Languages}{}%
  \ifthenelse{\equal{\@journal}{laws}}{Laws}{}%
  \ifthenelse{\equal{\@journal}{life}}{Life}{}%
  \ifthenelse{\equal{\@journal}{literature}}{Literature}{}%
  \ifthenelse{\equal{\@journal}{lubricants}}{Lubricants}{}%
  \ifthenelse{\equal{\@journal}{machines}}{Machines}{}%
  \ifthenelse{\equal{\@journal}{magnetochemistry}}{Magnetochemistry}{}%
  \ifthenelse{\equal{\@journal}{marinedrugs}}{Mar. Drugs}{}%
  \ifthenelse{\equal{\@journal}{materials}}{Materials}{}%
  \ifthenelse{\equal{\@journal}{mathematics}}{Mathematics}{}%
  \ifthenelse{\equal{\@journal}{mca}}{Math. Comput. Appl.}{}%
  \ifthenelse{\equal{\@journal}{medsci}}{Med. Sci.}{}% 
  \ifthenelse{\equal{\@journal}{medicines}}{Medicines}{}%
  \ifthenelse{\equal{\@journal}{membranes}}{Membranes}{}%
  \ifthenelse{\equal{\@journal}{metabolites}}{Metabolites}{}%
  \ifthenelse{\equal{\@journal}{metals}}{Metals}{}%
  \ifthenelse{\equal{\@journal}{microarrays}}{Microarrays}{}%
  \ifthenelse{\equal{\@journal}{micromachines}}{Micromachines}{}%
  \ifthenelse{\equal{\@journal}{microorganisms}}{Microorganisms}{}%
  \ifthenelse{\equal{\@journal}{minerals}}{Minerals}{}%
  \ifthenelse{\equal{\@journal}{molbank}}{Molbank}{}%
  \ifthenelse{\equal{\@journal}{molecules}}{Molecules}{}%
  \ifthenelse{\equal{\@journal}{mps}}{Methods and Protocols}{}%
  \ifthenelse{\equal{\@journal}{nanomaterials}}{Nanomaterials}{}%
  \ifthenelse{\equal{\@journal}{ncrna}}{Non-coding RNA}{}%
  \ifthenelse{\equal{\@journal}{neonatalscreening}}{Int. J. Neonatal Screen.}{}%
  \ifthenelse{\equal{\@journal}{nutrients}}{Nutrients}{}%
  \ifthenelse{\equal{\@journal}{particles}}{Particles}{}%
  \ifthenelse{\equal{\@journal}{pathogens}}{Pathogens}{}%
  \ifthenelse{\equal{\@journal}{pharmaceuticals}}{Pharmaceuticals}{}%
  \ifthenelse{\equal{\@journal}{pharmaceutics}}{Pharmaceutics}{}%
  \ifthenelse{\equal{\@journal}{pharmacy}}{Pharmacy}{}%
  \ifthenelse{\equal{\@journal}{philosophies}}{Philosophies}{}%
  \ifthenelse{\equal{\@journal}{photonics}}{Photonics}{}%
  \ifthenelse{\equal{\@journal}{plants}}{Plants}{}%
  \ifthenelse{\equal{\@journal}{polymers}}{Polymers}{}%
  \ifthenelse{\equal{\@journal}{processes}}{Processes}{}%
  \ifthenelse{\equal{\@journal}{proteomes}}{Proteomes}{}%
  \ifthenelse{\equal{\@journal}{publications}}{Publications}{}%
  \ifthenelse{\equal{\@journal}{recycling}}{Recycling}{}%
  \ifthenelse{\equal{\@journal}{religions}}{Religions}{}%
  \ifthenelse{\equal{\@journal}{remotesensing}}{Remote Sens.}{}%
  \ifthenelse{\equal{\@journal}{resources}}{Resources}{}%
  \ifthenelse{\equal{\@journal}{risks}}{Risks}{}%
  \ifthenelse{\equal{\@journal}{robotics}}{Robotics}{}%
  \ifthenelse{\equal{\@journal}{safety}}{Safety}{}%
  \ifthenelse{\equal{\@journal}{sensors}}{Sensors}{}%
  \ifthenelse{\equal{\@journal}{sexes}}{Sexes}{}%
  \ifthenelse{\equal{\@journal}{sinusitis}}{Sinusitis}{}%
  \ifthenelse{\equal{\@journal}{socsci}}{Soc. Sci.}{}%
  \ifthenelse{\equal{\@journal}{societies}}{Societies}{}%
  \ifthenelse{\equal{\@journal}{soils}}{Soils}{}%
  \ifthenelse{\equal{\@journal}{sports}}{Sports}{}%
  \ifthenelse{\equal{\@journal}{standards}}{Standards}{}%
  \ifthenelse{\equal{\@journal}{sustainability}}{Sustainability}{}%
  \ifthenelse{\equal{\@journal}{symmetry}}{Symmetry}{}%
  \ifthenelse{\equal{\@journal}{systems}}{Systems}{}%
  \ifthenelse{\equal{\@journal}{technologies}}{Technologies}{}%
  \ifthenelse{\equal{\@journal}{toxics}}{Toxics}{}%
  \ifthenelse{\equal{\@journal}{toxins}}{Toxins}{}%
  \ifthenelse{\equal{\@journal}{universe}}{Universe}{}%
  \ifthenelse{\equal{\@journal}{urbansci}}{Urban Sci.}{}%
  \ifthenelse{\equal{\@journal}{vaccines}}{Vaccines}{}%
  \ifthenelse{\equal{\@journal}{vetsci}}{Vet. Sci.}{}%
  \ifthenelse{\equal{\@journal}{viruses}}{Viruses}{}%
  \ifthenelse{\equal{\@journal}{visions}}{Visions}{}%
  \ifthenelse{\equal{\@journal}{water}}{Water}{}%
  }
    
%%%% Header and footer on first page
%% The plain page style needs to be redefined because with \maketitle in the article class, LaTeX applies the the plain page style automatically to the first page.
\ifthenelse{\equal{\@arttype}{Supfile}}{
\fancypagestyle{plain}{
\fancyhf{}
\fancyhead[R]{
\fontsize{8}{8}\selectfont S\thepage{} of S\pageref*{LastPage}}
\fancyhead[L]{\fontsize{8}{8}\selectfont
\ifthenelse{\equal{\@status}{submit}}{Version {\@ \today} submitted to {\em \journalname}}{{\em \journalname} {\bfseries \@pubyear}, {\em \@pubvolume}, \ifthenelse{\equal{\@articlenumber}{}}{\@firstpage --\pageref*{LastPage}}{\@articlenumber}; doi:\@doinum}
}}
}
{
\fancypagestyle{plain}{
\fancyhf{}
\fancyfoot[L]{
\fontsize{8}{8}\selectfont
  {\ifthenelse{\equal{\@status}{submit}}{\normalfont Submitted to {\em \journalname}, pages \thepage \ -- \color{black}{\pageref*{LastPage}}}
  {\normalfont {\em \journalname}\ {\bfseries {\@pubyear}}, {\em \@pubvolume}, \ifthenelse{\equal{\@articlenumber}{}}{\@firstpage --\pageref*{LastPage}}{\@articlenumber}; doi:\@doinum}}
}
 \fancyfoot[R]{\fontsize{8}{8}\selectfont{\changeurlcolor{black}\href{http://www.mdpi.com/journal/\@journal}{www.mdpi.com/journal/\@journal}}}
 \fancyhead{}
 \renewcommand{\headrulewidth}{0.0pt}%
}}

%%%% Maketitle
\renewcommand{\@maketitle}{
\ifthenelse{\equal{\@arttype}{Supfile}}{
\fontsize{18}{18}\selectfont
\noindent\textbf{Supplementary Materials: \@Title}
\par
\vspace{12pt}
\fontsize{10}{10}\selectfont
\noindent\boldmath\bfseries{\@Author}
}{
\begin{flushleft}
\vspace*{-1.75cm}
{%0
\ifthenelse{\equal{\@status}{submit}}{
\hfill \href{http://www.mdpi.com}{\includegraphics[height=1cm]{logo-mdpi}}\vspace{0.5cm}}
{\href{http://www.mdpi.com/journal/\@journal}{\includegraphics[height=1.2cm]{\@journal-logo}} \hfill \href{http://www.mdpi.com}{\includegraphics[height=1cm]{logo-mdpi}} \hspace{-6pt}}
\par
}%0
{%1
\vspace{14pt}
\fontsize{10}{10}\selectfont
\textit{\@arttype}
\par
}%1
{%2
\vspace{-1pt}
\fontsize{18}{18}\selectfont
\bfseries{\@Title}
\par
\vspace{15pt}
}%2
{%3
\fontsize{10}{10}\selectfont
\boldmath\bfseries{\@Author}
\par
\vspace{-4pt}
}%3

\end{flushleft}
}
}

\newcommand{\maketitlen}{
\begin{flushleft}
	\begin{spacing}{1.35}
{%4
\leftskip0.2cm
\fontsize{9}{9}\selectfont
\ifthenelse{\equal{\firstargument}{1}}{}{\hangafter=1\hangindent=1.7em \noindent}%
\@address
\par
}%4
{%5
\leftskip0.2cm
\fontsize{9}{9}\selectfont
\ifthenelse{\equal{\@authornum}{authors}}{
\ifthenelse{\equal{\@corres}{\@empty}}{}{\hangafter=1\hangindent=1.7em\noindent \textbf{*} \quad \@corres}
\par
}{}
}%5
{%6
\leftskip0.2cm
\fontsize{9}{9}\selectfont
\ifthenelse{\equal{\@firstnote}{\@empty}}{}{\hangafter=1\hangindent=1.7em\noindent  
$\dagger$ \quad \@firstnote}
\par
}%6
{%7
\leftskip0.2cm
\fontsize{9}{9}\selectfont
\ifthenelse{\equal{\@secondnote}{\@empty}}{}{\hangafter=1\hangindent=1.7em \noindent $\ddagger$ \quad \@secondnote}
\par
}%7
{%8
\leftskip0.2cm
\fontsize{9}{9}\selectfont
\ifthenelse{\equal{\@conference}{\@empty}}{}{\hangafter=1\hangindent=1.7em\noindent \textsection
 \quad This paper is an extended version of our paper published in\space \@conference.}
\par
}%8
{%9
\vspace{6pt}
\leftskip0.2cm
\fontsize{9}{9}\selectfont
\ifthenelse{\equal{\@externaleditor}{\@empty}}{}{\@externaleditor}
\par
}%9
{%10
\leftskip0.2cm
\fontsize{9}{9}\selectfont 
\ifthenelse{\equal{\@status}{submit}}{Version {\@ \today} submitted to \journalname; Typeset by \LaTeX{} using class file mdpi.cls}{\@history}
\par
\vspace{-4pt}
}%10
\end{spacing}
\end{flushleft}
}

%%%% Abstract, keywords, journal data, PACS, MSC, JEL
\newcommand{\abstractkeywords}{%
{%10
\begingroup
\vspace{-15pt}
\leftskip0.2cm 
\ifthenelse{\equal{\@simplesumm}{\@empty}}{}{
  \noindent\textbf{Simple Summary:\space}\@simplesumm
  \vspace{12pt}
  \par
}
\ifthenelse{\equal{\@abstract}{\@empty}}{}{
\noindent\textbf{Abstract:\space}\@abstract
\vspace{12pt}
\par
}
\endgroup
}%10
{%11
\begingroup
\leftskip0.2cm
\ifthenelse{\equal{\@keyword}{\@empty}}{}{
\noindent\textbf{Keywords:\space}\@keyword
\vspace{12pt}
\par
}
\endgroup
}%11
{% For journal Data:
\ifthenelse{\equal{\@dataset}{\@empty}}{}{
\begingroup
\leftskip0.2cm
\noindent\textbf{Data Set:\space}\@dataset
\vspace{12pt}
\par
\endgroup}
}%
{%For journal Data:
\ifthenelse{\equal{\@datasetlicense}{\@empty}}{}{
\begingroup
\leftskip0.2cm
\noindent\textbf{Data Set License:\space}\@datasetlicense
\vspace{12pt}
\par
\endgroup}
}%
{%12
\ifthenelse{\equal{\@pacs}{\@empty}}{}{
\begingroup
\leftskip0.2cm
\noindent\textbf{PACS:\space}\@pacs
\vspace{12pt}
\par
\endgroup}
}%12
{%13
\ifthenelse{\equal{\@msc}{\@empty}}{}{
\begingroup
\leftskip0.2cm
\noindent\textbf{MSC:\space}\@msc
\vspace{12pt}
\par
\endgroup}
}%13
{%14
\ifthenelse{\equal{\@jel}{\@empty}}{}{
\begingroup
\leftskip0.2cm
\noindent\textbf{JEL:\space}\@jel
\vspace{12pt}
\par
\endgroup}
}%14
\vspace{4pt}
\hrule
\vspace{12pt}
}

%%%% Print maketitle and abstractkeywords
\ifthenelse{\equal{\@arttype}{Supfile}}{
\AfterEndPreamble{
\maketitle
\let\maketitle\relax
\ifthenelse{\equal{\@status}{submit}}{\linenumbers}{}
}
}
{\AfterEndPreamble{
\maketitle
\let\maketitle\relax
\maketitlen
\let\maketitlen\relax
\ifthenelse{\equal{\@status}{submit}}{\linenumbers}{}
\abstractkeywords
}
}
\AtBeginDocument{
  \DeclareSymbolFont{AMSb}{U}{msb}{m}{n}
  \DeclareSymbolFontAlphabet{\mathbb}{AMSb}}

%%%% Section headings
\titleformat {\section} [block] {\raggedright \fontsize{10}{10}\selectfont\bfseries} {\thesection. \space} {0pt} {}
\titlespacing {\section} {0pt} {12pt} {6pt}

\titleformat {\subsection} [block] {\raggedright \fontsize{10}{10}\selectfont\itshape} {\thesubsection .\space} {0pt} {}
\titlespacing {\subsection} {0pt} {12pt} {6pt}

\titleformat {\subsubsection} [block] {\raggedright \fontsize{10}{10}\selectfont} {\thesubsubsection .\space} {0pt} {}
\titlespacing {\subsubsection} {0pt} {12pt} {6pt}

\titleformat {\paragraph} [block] {\raggedright \fontsize{10}{10}\selectfont} {} {0pt} {}
\titlespacing {\paragraph} {0pt} {12pt} {6pt}

%%%% Special section title style for back matter
\newcommand{\supplementary}[1]{
\par\vspace{6pt}\noindent{\fontsize{9}{9}\selectfont\textbf{Supplementary Materials:} {#1}\par}}

\newcommand{\acknowledgments}[1]{
\vspace{6pt}\noindent{\fontsize{9}{9}\selectfont\textbf{Acknowledgments:} {#1}\par}}

\newcommand{\authorcontributions}[1]{%
\vspace{6pt}\noindent{\fontsize{9}{9}\selectfont\textbf{Author Contributions:} {#1}\par}}

\newcommand{\conflictofinterests}[1]{%
\vspace{6pt}\noindent{\fontsize{9}{9}\selectfont\textbf{Conflicts of Interest:} {#1}\par}}

\newcommand{\sampleavailability}[1]{%
\vspace{6pt}\noindent{\fontsize{9}{9}\selectfont\textbf{Sample Availability:} {#1}\par}}

\newcommand{\abbreviations}[1]{%
\vspace{12pt}\noindent{\selectfont\textbf{Abbreviations}\par\vspace{6pt}\noindent {#1}\par}}

%%%%% Defines the appendix
\renewcommand{\appendix}{
\setcounter{section}{0}%
\setcounter{subsection}{0}%
\setcounter{subsubsection}{0}%
\gdef\thesection{\@Alph\c@section}%
\gdef\thesubsection{\@Alph\c@section.\@arabic\c@subsection}%
\titleformat {\section} [block] {\raggedright \fontsize{10}{10}\selectfont\bfseries} {\appendixname~\thesection.~} {0pt} {}
\titlespacing {\section} {0pt} {12pt} {6pt}
\titleformat {\subsection} [block] {\raggedright \fontsize{10}{10}\selectfont\itshape} {\appendixname~\thesubsection.~} {0pt} {}
\titlespacing {\subsection} {0pt} {12pt} {6pt}
\titleformat {\subsubsection} [block] {\raggedright \fontsize{10}{10}\selectfont} {\appendixname~\thesubsubsection.~} {0pt} {}
\titlespacing {\subsubsection} {0pt} {12pt} {6pt}
\gdef\theHsection{\@Alph\c@section.}% for hyperref
\gdef\theHsubsection{\@Alph\c@section.\@arabic\c@subsection}% for hyperref
\csname appendixmore\endcsname
}

%%%% Layout
\RequirePackage[left=2.7cm,right=2.7cm,top=1.8cm,bottom=1.5cm, includehead,includefoot]{geometry}
\linespread{1.13} 
\setlength{\parindent}{0.75cm}

%%%% Figures and tables
\newfloat{boxenv}{H}{lop}
\floatname{boxenv}{Box}
\floatstyle{plain}
\restylefloat{figure}
\floatstyle{plaintop}
\restylefloat{table}
\floatstyle{plaintop}
\restylefloat{boxenv}
\floatstyle{plain}
\restylefloat{scheme}
\RequirePackage{caption} 
\captionsetup[figure]{position=bottom, margin=0.75cm, labelfont={bf, small, stretch=1.17}, labelsep=period, textfont={small, stretch=1.17}, aboveskip=6pt, belowskip=12pt, justification=justified}
\captionsetup[scheme]{position=bottom, margin=0.75cm, labelfont={bf, small, stretch=1.17}, labelsep=period, textfont={small, stretch=1.17}, aboveskip=6pt, belowskip=12pt, justification=justified}
\captionsetup[table]{position=top, margin=0.75cm, labelfont={bf, small, stretch=1.17}, labelsep=period, textfont={small, stretch=1.17}, aboveskip=6pt, belowskip=12pt, justification=justified}
\captionsetup[boxenv]{position=top, margin=0.75cm, labelfont={bf, small, stretch=1.17}, labelsep=period, textfont={small, stretch=1.17}, aboveskip=6pt, belowskip=12pt, justification=justified}

%% For table footnotes
\newsavebox{\@justcentbox}
\newcommand{\justifyorcenter}[1]{
\sbox \@justcentbox{#1}
\ifdim \wd \@justcentbox >\hsize #1
\else \centerline{#1} \fi
}

%%%% Bullet lists
\newlength{\wideitemsep}
\setlength{\wideitemsep}{.5\itemsep}
\addtolength{\wideitemsep}{-7pt}
\let\olditem\item
\renewcommand{\item}{\setlength{\itemsep}{\wideitemsep}\olditem}

%%%% Supplementary file
\ifthenelse{\equal{\@arttype}{Supfile}}
{%
  \renewcommand{\thefigure}{S\arabic{figure}}%
  \renewcommand{\thetable}{S\arabic{table}}%
}{}%

%%%% Header and footer (all pages except the first)
\renewcommand\headrule{} %% set line (from fancyhdr) in header to nothing
\pagestyle{fancy}
\lhead{
\fontsize{8}{8}\selectfont
\ifthenelse{\equal{\@status}{submit}}{Version {\@ \today} submitted to {\em \journalname}}{
\ifthenelse{\equal{\@arttype}{Supfile}}{{\em \journalname} {\bfseries \@pubyear}, {\em \@pubvolume}, \ifthenelse{\equal{\@articlenumber}{}}{\@firstpage --\pageref*{LastPage}}{\@articlenumber}; doi:\@doinum}{
{\em \journalname} {\bfseries \@pubyear}, {\em \@pubvolume}, \ifthenelse{\equal{\@articlenumber}{}}{\@firstpage --\pageref*{LastPage}}{\@articlenumber}}}
}
\ifthenelse{\equal{\@journal}{recycling} \OR \equal{\@journal}{philosophies} \OR \equal{\@journal}{sinusitis} \OR \equal{\@journal}{fishes} \OR \equal{\@journal}{jfmk}}{
\rhead{
\ifthenelse{\equal{\@arttype}{Supfile}}{\fontsize{8}{8}\selectfont S\thepage{}}{
\ifthenelse{\equal{\@status}{submit}}{\fontsize{8}{8}\selectfont\thepage}{
\fontsize{8}{8}\selectfont\thepage{}}}
}
}
{
\rhead{
\ifthenelse{\equal{\@arttype}{Supfile}}{
\fontsize{8}{8}\selectfont S\thepage{} of S\pageref*{LastPage}
}
{
\ifthenelse{\equal{\@status}{submit}}{\fontsize{8}{8}\selectfont\thepage \ of \pageref*{LastPage}}{
\fontsize{8}{8}\selectfont\thepage{} of \pageref*{LastPage}}
}}
}
\cfoot{}

%%%% Bibliography
\renewcommand\@biblabel[1]{#1.\hfill}
\def\thebibliography#1{
\linespread{1.44} 
\vspace{12pt}
\fontsize{10}{10}\selectfont
\noindent{{\bfseries \bibname}}
\vspace{6pt}
\fontsize{9}{9}\selectfont
\list{{\arabic{enumi}}}{\def\makelabel##1{\hss{##1}}
\topsep=0\p@
\parsep=5\p@
\partopsep=0\p@
\itemsep=0\p@
\labelsep=1.5mm
\itemindent=0\p@
\settowidth\labelwidth{\footnotesize[#1]}%
\leftmargin\labelwidth
\advance\leftmargin\labelsep
%\advance\leftmargin -\itemindent
\usecounter{enumi}}
%\def\newblock{\ }
%\sloppy\clubpenalty4000\widowpenalty4000
%\sfcode`\.=1000\relax
}
\let\endthebibliography=\endlist

%%%% Copyright info
\newcommand{\cright}{%
\ifthenelse{\equal{\@arttype}{Supfile}}{}{%
\vspace{12pt}
\noindent
\linespread{1.44} 
\fontsize{9}{9}\selectfont 
\ifthenelse{\equal{\@status}{submit}}{\noindent \copyright{} {\@ \the\year} by the
\@authornum. Submitted to {\em \journalname} for possible open
access publication under the terms and conditions of the
Creative Commons Attribution license
(\changeurlcolor{black}\href{http://creativecommons.org/licenses/by/4.0/.}{http://creativecommons.org/licenses/by/4.0/}) }{
\begin{minipage}{.2\textwidth}
\hspace{-1.2mm}\vspace{2mm}\href{http://creativecommons.org/}{\includegraphics[height=1cm]{logo-ccby}}
\end{minipage}% This must go next to `\end{minipage}`
\begin{minipage}{.8\textwidth}
\copyright
\ {\@copyrightyear} by
the \@authornum; licensee MDPI, Basel, Switzerland. This article is an
open access article distributed under the terms and conditions
of the Creative Commons by Attribution (CC-BY) license
(\changeurlcolor{black}\href{http://creativecommons.org/licenses/by/4.0/}{http://creativecommons.org/licenses/by/4.0/}).
\end{minipage}
}
}
}

\endinput