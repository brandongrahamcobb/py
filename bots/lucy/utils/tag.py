from typing import Optional, List, Dict

class TagManager:
    def __init__(self, db_pool):
        self.pool = db_pool

    async def add_tag(
        self,
        name: str,
        location_id: int,
        owner_id: int,
        content: Optional[str] = None,
        attachment_url: Optional[str] = None,
    ):
        query = """
            INSERT INTO tags (name, location_id, content, attachment_url, owner_id)
            VALUES ($1, $2, $3, $4, $5)
        """
        async with self.pool.acquire() as conn:
            await conn.execute(query, name, location_id, content, attachment_url, owner_id)

    async def get_tag(self, location_id: int, name: str) -> Dict:
        query = """SELECT * FROM tags WHERE location_id = $1 AND LOWER(name) = $2"""
        async with self.pool.acquire() as conn:
            tag = await conn.fetchrow(query, location_id, name.lower())
            if tag:
                return dict(tag)
            raise RuntimeError(f'Tag "{name}" not found.')

    async def update_tag(
        self,
        name: str,
        location_id: int,
        owner_id: int,
        content: Optional[str] = None,
        attachment_url: Optional[str] = None,
    ) -> int:
        query = """
            UPDATE tags
            SET content = $1, attachment_url = $2
            WHERE name = $3 AND location_id = $4 AND owner_id = $5
        """
        async with self.pool.acquire() as conn:
            result = await conn.execute(query, content, attachment_url, name, location_id, owner_id)
            return int(result.split(" ")[-1])  # Extract the number of rows affected

    async def delete_tag(self, name: str, location_id: int, owner_id: int) -> int:
        query = """DELETE FROM tags WHERE name = $1 AND location_id = $2 AND owner_id = $3"""
        async with self.pool.acquire() as conn:
            result = await conn.execute(query, name, location_id, owner_id)
            return int(result.split(" ")[-1])  # Extract the number of rows affected

    async def list_tags(self, location_id: int, owner_id: int) -> List[Dict]:
        query = """SELECT name, content, attachment_url FROM tags WHERE location_id = $1 AND owner_id = $2"""
        async with self.pool.acquire() as conn:
            tags = await conn.fetch(query, location_id, owner_id)
        return [dict(tag) for tag in tags]
