.PHONY: up down logs migrate etl-npa dump restore

up:
		docker compose up -d db adminer

down:
		docker compose down

logs:
		docker compose logs -f db

migrate:
		python -m bionexus.cli upgrade head

etl-npa:
		python -m bionexus.cli load-npatlas --path ./data/npatlas.csv --format auto --dedupe-on inchikey

dump:
		python -m bionexus.cli dump-db --out dumps/bionexus_`date +%Y%m%d`.dump

restore:
		python -m bionexus.cli restore-db --dump $(dumpfile)
