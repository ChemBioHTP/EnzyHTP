#! python3
# -*- encoding: utf-8 -*-
'''
This file is used to define the table structure/model of the database 
that records execution information for each WorkUnit and WorkFlow instance.

@File    :   database.py
@Created :   2024/03/17 14:40
@Author  :   Zhong, Yinjie
@Version :   1.0
@Contact :   yinjie.zhong@vanderbilt.edu
'''

# Here put the import lib.
from sqlalchemy import Column, Integer, String, DateTime, create_engine
from sqlalchemy.orm import sessionmaker, Session, declarative_base
from datetime import datetime
from os import path, remove

# Declare base object and connections.
Base = declarative_base()

class SqliteSqlalchemy(object):
    """
    A class for managing the SQLite database connection and session using SQLAlchemy.

    Attributes:
        session (Session): An instance of SQLAlchemy session for database operations.

    Methods:
        __init__: Constructor for the class, establishes the database connection and session.
    """

    def __init__(self, db_filepath: str, overwrite: bool = False):
        """
        Initializes the database connection and session.

        Args:
            db_filepath (str): The file path for the SQLite database.
            overwrite: If overwrite the database file when it exists. Default False.
        """
        if (overwrite and path.isfile(db_filepath)):
            remove(path=db_filepath)
        engine = create_engine(f"sqlite:///{db_filepath}", echo=False)  # Create connection engine.
        Base.metadata.create_all(engine, checkfirst=True)   # Create tables.
        self.session = sessionmaker(bind=engine)()  # Create session.

# Define models.
class ExecutionEntity(Base):
    """
    A SQLAlchemy ORM model for storing and managing execution entities in the database.

    Attributes:
        identifier (str): A unique identifier for the execution entity.
        status (int): The execution status of the entity.
        updated_time (datetime): The timestamp of the last update to the entity.

    Methods:
        __init__: Constructor for creating an ExecutionEntity instance.
        insert_or_update: Inserts a new entity or updates an existing one in the database.
        get_by_identifier: Retrieves an entity by its unique identifier from the database.
        get_all: Retrieves all entities from the database.
    """
    __tablename__ = "execution_entities"

    identifier = Column(String(255), primary_key=True, unique=True)
    status = Column(Integer, nullable=False)
    updated_time = Column(DateTime, nullable=False)

    def __repr__(self) -> str:
        cls = self.__class__.__name__
        return f"{cls}(identifier={self.identifier}, status={self.status}, updated_time={self.updated_time})"
    
    def __init__(self, identifier: str, status: int):
        """
        Initializes a new instance of ExecutionEntity.

        Args:
            identifier (str): The unique identifier for the entity.
            status (int): The execution status of the entity.
        """
        self.identifier = identifier
        self.status = status
        self.updated_time = datetime.now()
        return
    
    def insert_or_update(self, session: Session):
        """
        Inserts a new entity into the database or updates an existing one based on the identifier.

        Args:
            session (Session): The SQLAlchemy session for database operations.
        """
        existing_entity = session.query(ExecutionEntity).filter_by(identifier=self.identifier).first()
        if existing_entity:
            existing_entity.status = self.status
            existing_entity.updated_time = datetime.now()
        else:
            session.add(self)
        session.commit()

    @staticmethod
    def get_by_identifier(session: Session, identifier: str):
        """
        Retrieves an entity from the database based on its unique identifier.

        Args:
            session (Session): The SQLAlchemy session for database operations.
            identifier (str): The unique identifier of the entity to retrieve.

        Returns:
            ExecutionEntity: The retrieved entity, if found; otherwise, None.
        """
        return session.query(ExecutionEntity).filter_by(identifier=identifier).first()

    @staticmethod
    def get_all(session):
        """
        Retrieves all entities from the database.

        Args:
            session (Session): The SQLAlchemy session for database operations.

        Returns:
            list: A list of all ExecutionEntity instances in the database.
        """
        return session.query(ExecutionEntity).all()
    